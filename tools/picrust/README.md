## PICRUSt Galaxy Wrapper

This page is for anyone interested in installing [PICRUSt](http://picrust.github.io/picrust/) in a [galaxy instance](https://galaxyproject.org/).The below descriptions should clarify the most confusing part of the installation.

The most important steps are that you need to go to the Galaxy IUC tool shed and install the 4 wrapper scripts there along with PICRUSt itself. If you install PICRUSt through this interface the precalculated files will not also be installed, you will need to download them separately (see below).

Once these scripts are installed sections like the one below should be added to *galaxy/config/shed_tool_data_table_conf.xml* automatically. They will point to different *picrust_precalculated.loc* files; you only need to edit one of them since all of these files will be checked.

```
<tables>
    <table name="picrust_precalculated" comment_char="#">
        <columns>name, value</columns>
        <file path="tool-data/picrust_precalculated.loc" />
    </table>
</tables>
```

You should download the precalculated files from http://kronos.pharmacology.dal.ca/public_files/picrust/picrust_precalculated_v1.1.1/13_5/.

Most users only use the below two files (however you can also download the precalculated COG and RFAM abundances):
```
16S_13_5_precalculated.tab.gz
ko_13_5_precalculated.tab.gz
```

Once you download these files we just need to point to them in a *picrust_precalculated.loc* file.

An example file is shown below. Note that tabs should separate the two columns and that in the below command "16S 13_5" (separated by one space) is one example dbname.

```
<dbname>       <file_base>
16S 13_5    /path/to/16S_13_5_precalculated.tab.gz
ko 13_5    /path/to/ko_13_5_precalculated.tab.gz
```
