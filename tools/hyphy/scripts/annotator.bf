LoadFunctionLibrary ("libv3/tasks/trees.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");
LoadFunctionLibrary ("libv3/IOFunctions.bf");

tags = {
    "^Node" : "Internal"
};




tree = trees.LoadAnnotatedTopology (TRUE);
ts = tree[^"terms.trees.newick_with_lengths"];
root_node = io.PromptUserForString ("Root on this node");



Topology T = ts;
ts = RerootTree (T, root_node);
ACCEPT_ROOTED_TREES = FALSE;
Topology T = ts;
ACCEPT_ROOTED_TREES = TRUE;
Topology TR = ts;


NORMALIZE_SEQUENCE_NAMES = FALSE;
SetDialogPrompt ("File with the sequences to label as the in-clade");
DataSet query = ReadDataFile (PROMPT_FOR_FILE);
GetString (seqNames,query,-1);

label = io.PromptUserForString ("Use this label");


for (s;in;seqNames) {
    tags [s && 6] = label;
}

reg_exp = Rows (tags);


node_labels = {};
for (_regexp_, _leaves_; in; regexp.PartitionByRegularExpressions (BranchName (T,-1), reg_exp)) {
    tag = tags[_regexp_];
    
    if (tag != "Internal") {    
    
        if (Abs (tag) == 0) {
            tag = "Reference";
        }
    
        for (l; in; _leaves_) {
            node_labels[l] = tag;
        }
    }
}

node_labelsR = {};
for (_regexp_, _leaves_; in; regexp.PartitionByRegularExpressions (BranchName (TR,-1), reg_exp)) {
    tag = tags[_regexp_];
    
    if (tag != "Internal") {    
    
        if (Abs (tag) == 0) {
            tag = "Reference";
        }
    
        for (l; in; _leaves_) {
            node_labelsR[l] = tag;
        }
    }
}


leaf_labels = node_labels;
node_labelsF = node_labels;
node_labels * ((trees.ParsimonyLabel ("T", node_labels))["labels"]);
node_labelsR * ((trees.ParsimonyLabel ("T", node_labelsR))["labels"]);
node_labelsF * ((trees.ParsimonyLabel ("T", node_labelsF))["labels"]);

output_to = io.PromptUserForString ("Write output to this prefix");

fprintf (output_to + "labels.json", CLEAR_FILE, leaf_labels);
fprintf (output_to + "int.nwk", CLEAR_FILE, tree.Annotate ("T", "relabel_and_annotate", "{}", FALSE));
fprintf (output_to + "clade.nwk", CLEAR_FILE, tree.Annotate ("TR", "relabel_and_annotate_full", "{}", FALSE));
fprintf (output_to + "full.nwk", CLEAR_FILE, tree.Annotate ("T", "relabel_and_annotate_full", "{}", FALSE));

function relabel_and_annotate (node_name) {
    _label = "";
    if (node_labels / node_name && leaf_labels / node_name == FALSE) {
        _label = "{" + node_labels[node_name] + "}";
    }
    return node_name + _label;
}

function relabel_and_annotate_full (node_name) {
    _label = "";
    if (node_labelsR / node_name) {
        _label = "{" + node_labelsR[node_name] + "}";
    }
    return node_name + _label;
}
