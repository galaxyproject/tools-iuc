#!/usr/bin/env python
import argparse
import copy
import logging
import sys

from BCBio import GFF
from Bio.SeqFeature import FeatureLocation

logging.basicConfig(level=logging.INFO)
log = logging.getLogger(__name__)

__author__ = "Eric Rasche"
__version__ = "0.4.0"
__maintainer__ = "Eric Rasche"
__email__ = "esr@tamu.edu"


def feature_lambda(feature_list, test, test_kwargs, subfeatures=True):
    """Recursively search through features, testing each with a test function, yielding matches.

    GFF3 is a hierachical data structure, so we need to be able to recursively
    search through features. E.g. if you're looking for a feature with
    ID='bob.42', you can't just do a simple list comprehension with a test
    case. You don't know how deeply burried bob.42 will be in the feature tree. This is where feature_lambda steps in.

    :type feature_list: list
    :param feature_list: an iterable of features

    :type test: function reference
    :param test: a closure with the method signature (feature, **kwargs) where
                 the kwargs are those passed in the next argument. This
                 function should return True or False, True if the feature is
                 to be yielded as part of the main feature_lambda function, or
                 False if it is to be ignored. This function CAN mutate the
                 features passed to it (think "apply").

    :type test_kwargs: dictionary
    :param test_kwargs: kwargs to pass to your closure when it is called.

    :type subfeatures: boolean
    :param subfeatures: when a feature is matched, should just that feature be
                        yielded to the caller, or should the entire sub_feature
                        tree for that feature be included? subfeatures=True is
                        useful in cases such as searching for a gene feature,
                        and wanting to know what RBS/Shine_Dalgarno_sequences
                        are in the sub_feature tree (which can be accomplished
                        with two feature_lambda calls). subfeatures=False is
                        useful in cases when you want to process (and possibly
                        return) the entire feature tree, such as applying a
                        qualifier to every single feature.

    :rtype: yielded list
    :return: Yields a list of matching features.
    """
    # Either the top level set of [features] or the subfeature attribute
    for feature in feature_list:
        if test(feature, **test_kwargs):
            if not subfeatures:
                feature_copy = copy.deepcopy(feature)
                feature_copy.sub_features = []
                yield feature_copy
            else:
                yield feature

        if hasattr(feature, 'sub_features'):
            for x in feature_lambda(feature.sub_features, test, test_kwargs, subfeatures=subfeatures):
                yield x


def feature_test_qual_value(feature, **kwargs):
    """Test qualifier values.

    For every feature, check that at least one value in
    feature.quailfiers(kwargs['qualifier']) is in kwargs['attribute_list']
    """
    for attribute_value in feature.qualifiers.get(kwargs['qualifier'], []):
        if attribute_value in kwargs['attribute_list']:
            return True
    return False


def __get_features(child, interpro=False):
    child_features = {}
    for rec in GFF.parse(child):
        # Only top level
        for feature in rec.features:
            # Get the record id as parent_feature_id (since this is how it will be during remapping)
            parent_feature_id = rec.id
            # If it's an interpro specific gff3 file
            if interpro:
                # Then we ignore polypeptide features as they're useless
                if feature.type == 'polypeptide':
                    continue
                # If there's an underscore, we strip up to that underscore?
                # I do not know the rationale for this, removing.
                # if '_' in parent_feature_id:
                    # parent_feature_id = parent_feature_id[parent_feature_id.index('_') + 1:]

            try:
                child_features[parent_feature_id].append(feature)
            except KeyError:
                child_features[parent_feature_id] = [feature]
            # Keep a list of feature objects keyed by parent record id
    return child_features


def __update_feature_location(feature, parent, protein2dna):
    start = feature.location.start
    end = feature.location.end
    if protein2dna:
        start *= 3
        end *= 3

    if parent.location.strand >= 0:
        ns = parent.location.start + start
        ne = parent.location.start + end
        st = +1
    else:
        ns = parent.location.end - end
        ne = parent.location.end - start
        st = -1

    # Don't let start/stops be less than zero. It's technically valid for them
    # to be (at least in the model I'm working with) but it causes numerous
    # issues.
    #
    # Instead, we'll replace with %3 to try and keep it in the same reading
    # frame that it should be in.
    if ns < 0:
        ns %= 3
    if ne < 0:
        ne %= 3

    feature.location = FeatureLocation(ns, ne, strand=st)

    if hasattr(feature, 'sub_features'):
        for subfeature in feature.sub_features:
            __update_feature_location(subfeature, parent, protein2dna)


def rebase(parent, child, interpro=False, protein2dna=False, map_by='ID'):
    # get all of the features we will be re-mapping in a dictionary, keyed by parent feature ID
    child_features = __get_features(child, interpro=interpro)

    for rec in GFF.parse(parent):
        replacement_features = []
        for feature in feature_lambda(
                rec.features,
                # Filter features in the parent genome by those that are
                # "interesting", i.e. have results in child_features array.
                # Probably an unnecessary optimisation.
                feature_test_qual_value,
                {
                    'qualifier': map_by,
                    'attribute_list': child_features.keys(),
                },
                subfeatures=False):

            # Features which will be re-mapped
            to_remap = child_features[feature.id]
            # TODO: update starts
            fixed_features = []
            for x in to_remap:
                # Then update the location of the actual feature
                __update_feature_location(x, feature, protein2dna)

                if interpro:
                    for y in ('status', 'Target'):
                        try:
                            del x.qualifiers[y]
                        except Exception:
                            pass

                fixed_features.append(x)
            replacement_features.extend(fixed_features)
        # We do this so we don't include the original set of features that we
        # were rebasing against in our result.
        rec.features = replacement_features
        rec.annotations = {}
        GFF.write([rec], sys.stdout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('parent', type=argparse.FileType('r'), help='Parent GFF3 annotations')
    parser.add_argument('child', type=argparse.FileType('r'), help='Child GFF3 annotations to rebase against parent')
    parser.add_argument('--interpro', action='store_true',
                        help='Interpro specific modifications')
    parser.add_argument('--protein2dna', action='store_true',
                        help='Map protein translated results to original DNA data')
    parser.add_argument('--map_by', help='Map by key', default='ID')
    args = parser.parse_args()
    rebase(**vars(args))
