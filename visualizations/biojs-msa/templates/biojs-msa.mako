<%
	# based on scatterplot

    default_title = "BioJS MSA of '" + hda.name + "'"
    info = hda.name
    if hda.info:
        info += ' : ' + hda.info

    # Use root for resource loading.
    root = h.url_for( '/' )
%>
## ----------------------------------------------------------------------------

<!DOCTYPE HTML>
<html>
<head>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8" />
<title>${title or default_title} | ${visualization_display_name}</title>
${h.javascript_link( root + 'plugins/visualizations/biojs-msa/static/msa.min.js' )}

</head>

## ----------------------------------------------------------------------------
<body>
<div>
    <h2>${title or default_title}</h2>
    <p>${info}</p>
</div>

<div id="msa_menubar"></div>
<div id="msa_menu"></div>

<script type="text/javascript">
    var config  = ${h.dumps( config )};
    var hdaJson = ${h.dumps( trans.security.encode_dict_ids( hda.to_dict() ), indent=2 )};

    var title   = "${title or default_title}";
    var url = "${root}" + "/api/datasets/"+config.dataset_id+"?data_type=raw_data&provider=base";
	var xhr = require("nets");
	xhr(url, function(err, response,text){
		var data = JSON.parse(text).data;

        if(hdaJson.data_type == "galaxy.datatypes.sequence.Fasta"){
		  var seqs = require("biojs-io-fasta").parse.parse(data);
        }else{
          // it could be clustal
          var seqs = require("biojs-io-clustal").parse(data); 
        }

        var msa = require("biojs-vis-msa");
        
        var opts = {};
        
        opts.el = document.getElementById('msa_menu');
        opts.seqs = seqs;
        
        //opts.zoomer = { textVisible: false};
        opts.vis = {overviewbox: true};
        //opts.columns = {hidden: [1,2,3]};
        var m = new msa.msa(opts);
        
        // the menu is independent to the MSA container
        var menuOpts = {};
        menuOpts.el = document.getElementById("msa_menubar");
        menuOpts.msa = m;
        var defMenu = new msa.menu.defaultmenu(menuOpts);
        m.addView("menu", defMenu);

        m.render();
	});
</script>

</body>
</html>
