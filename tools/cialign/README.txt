Planemo Test: 
planemo test cialign.xml --galaxy_root /Users/ahmad777/project/galaxy --biocontainers
planemo test cialign.xml --galaxy_root /Users/ahmad777/project/galaxy  --conda_auto_init --conda_auto_install --update_test_data

--update_test_data              Update test-data directory with job outputs
--conda_auto_install / --no_conda_auto_install
                                Conda dependency resolution for Galaxy will
                                attempt to install requested but missing
                                packages

Order of parameter attributes:

    name
    argument
    type
    format
    min | truevalue
    max | falsevalue
    value | checked
    optional
    label
    help


Steps to include a parameter:
    1. help
    2. command
    3. paramater
    4. test
    

Archieve:
<param argument="--all" type="boolean" truevalue="--all" falsevalue="" label="Use all available functions, with default parameters unless others are specified."/>

<!-- Test 2 default value  (Basic) -->
        <test expect_num_outputs="3">
            <param name="input" value="example1.fasta"/>
        </test>
        <!-- Test 3 config file -->
        <test expect_num_outputs="3">
            <param name="input" value="example1.fasta"/>
            <param name="path_to_config" value="ini_template.ini"/>
            <output name="cleaned" file="CIAlign_cleaned.fasta"/>
            <output name="removed" file="CIAlign_removed.txt"/>
            <output name="log" file="CIAlign_log.txt"/>
        </test>
        <!-- Test 4 all options -->
        <test expect_num_outputs="3">
            <param name="input" value="example1.fasta"/>
            <param name="prefix" value="all"/>
            <param name="all_options" value="y"/>
            <output name="cleaned" file="all_cleaned.fasta"/>
            <output name="removed" file="all_removed.txt"/>
            <output name="log" file="all_log.txt"/>
        </test>
        <!-- Test 5 all clean options -->
        <test expect_num_outputs="3">
            <param name="input" value="example1.fasta"/>
            <param name="prefix" value="clean"/>
            <param name="clean" value="y"/>
            <output name="cleaned" file="clean_cleaned.fasta"/>
            <output name="removed" file="clean_removed.txt"/>
            <output name="log" file="clean_log.txt"/>
        </test>
        <!-- Test 6 all visualize and silent options -->
        <test expect_num_outputs="3">
            <param name="input" value="example1.fasta"/>
            <param name="prefix" value="visualise"/>
            <param name="clean" value="n"/>
            <param name="visualise" value="y"/>
            <param name="silent" value="y"/>
            <output name="cleaned" file="visualise_cleaned.fasta"/>
            <output name="removed" file="visualise_removed.txt"/>
            <output name="log" file="visualise_log.txt"/>
        </test>