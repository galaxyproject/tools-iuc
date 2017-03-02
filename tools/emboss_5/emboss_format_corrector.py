# EMBOSS format corrector
import operator


# Properly set file formats before job run
def exec_before_job( app, inp_data=None, out_data=None, tool=None, param_dict=None ):
    # why isn't items an ordered list?
    items = out_data.items()
    items = sorted(items, key=operator.itemgetter(0))

    # normal filetype correction
    data_count = 1
    for name, data in items:
        outputType = param_dict.get( 'out_format' + str(data_count), None )
        if outputType is not None:
            if outputType == 'ncbi':
                outputType = "fasta"
            elif outputType == 'excel':
                outputType = "tabular"
            elif outputType == 'text':
                outputType = "txt"
            data = app.datatypes_registry.change_datatype(data, outputType)
            app.model.context.add( data )
            app.model.context.flush()
        data_count += 1

    # html filetype correction
    data_count = 1
    for name, data in items:
        wants_plot = param_dict.get( 'html_out' + str(data_count), None )
        ext = "html"
        if wants_plot == "yes":
            data = app.datatypes_registry.change_datatype(data, ext)
            app.model.context.add( data )
            app.model.context.flush()
        data_count += 1

    # png file correction
    data_count = 1
    for name, data in items:
        wants_plot = param_dict.get( 'plot' + str(data_count), None )
        ext = "png"
        if wants_plot == "yes":
            data = app.datatypes_registry.change_datatype(data, ext)
            app.model.context.add( data )
            app.model.context.flush()
        data_count += 1
