def validate_input( trans, error_map, param_values, page_param_map ):
    """
        Validates the user input, before execution.
    """
    factors = param_values['rep_factorName']
    factor_name_list = []
    factor_duplication = False
    level_duplication = False

    for factor in factors:
        # factor names should be unique
        fn = factor['factorName']
        if fn in factor_name_list:
            factor_duplication = True
            break
        factor_name_list.append( fn )

        level_name_list = list()

        for level in ['factorLevel1', 'factorLevel2']:
            # level names under one factor should be unique
            fl = factor[level]
            if fl in level_name_list:
                level_duplication = True
            level_name_list.append( fl )

        if level_duplication:
            error_map['rep_factorName'] = [ dict() for t in factors ]
            for i in range( len( factors ) ):
                error_map['rep_factorName'][i]['FactorLevel1'] = [ {'factorLevel': 'Factor levels for each factor need to be unique'} for t in [factor['factorLevel1'], factor['factorLevel2']] ]
            break

    if factor_duplication:
        error_map['rep_factorName'] = [ dict() for t in factors ]
        for i in range( len( factors ) ):
            error_map['rep_factorName'][i]['factorName'] = 'Factor names need to be unique.'
