#!/usr/bin/env python
import sys
import logging
logging.basicConfig(level=logging.INFO)
import argparse
from BCBio import GFF
from Bio.SeqFeature import FeatureLocation
log = logging.getLogger(__name__)

__author__ = "Eric Rasche"
__version__ = "0.4.0"
__maintainer__ = "Eric Rasche"
__email__ = "esr@tamu.edu"

def __get_features(child, interpro=False):
    child_features = {}
    for rec in GFF.parse(child):
        for feature in rec.features:
            parent_feature_id = rec.id
            if interpro:
                parent_feature_id = parent_feature_id[parent_feature_id.index('_') + 1:]
                if feature.type == 'polypeptide':
                    continue

            try:
                child_features[parent_feature_id].append(feature)
            except KeyError:
                child_features[parent_feature_id] = [feature]
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


def rebase(parent, child, interpro=False, protein2dna=False):
    child_features = __get_features(child, interpro=interpro)

    for rec in GFF.parse(parent):
        # TODO, replace with recursion in case it's matched against a
        # non-parent feature. We're cheating a bit here right now...
        replacement_features = []
        for feature in rec.features:
            if feature.id in child_features:
                new_subfeatures = child_features[feature.id]
                # TODO: update starts
                fixed_subfeatures = []
                for x in new_subfeatures:
                    # Then update the location of the actual feature
                    __update_feature_location(x, feature, protein2dna)

                    if interpro:
                        for y in ('status', 'Target'):
                            try:
                                del x.qualifiers[y]
                            except:
                                pass

                    fixed_subfeatures.append(x)
                replacement_features.extend(fixed_subfeatures)
        # We do this so we don't include the original set of features that we
        # were rebasing against in our result.
        rec.features = replacement_features
        GFF.write([rec], sys.stdout)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='rebase gff3 features against parent locations', epilog="")
    parser.add_argument('parent', type=file, help='Parent GFF3 annotations')
    parser.add_argument('child', help='Child GFF3 annotations to rebase against parent')
    parser.add_argument('--interpro', action='store_true',
                        help='Interpro specific modifications')
    parser.add_argument('--protein2dna', action='store_true',
                        help='Map protein translated results to original DNA data')
    args = parser.parse_args()
    rebase(**vars(args))
