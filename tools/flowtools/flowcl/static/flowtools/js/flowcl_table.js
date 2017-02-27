// Copyright (c) 2016 Northrop Grumman.
// All rights reserved.

var flowcl_table = './flowCL_table.txt';

var displayFlowCLTable = function(){
  d3.tsv(flowcl_table, function(error, data){
    if (error){
      alert("Problem retrieving data");
      return;
    }
    var fclHTML = '<table id="flowcltable" class="cltable display compact" cellspacing="0" width="100%"/>';
    $('#flowcl').html(fclHTML);

    var fclTableData = $.extend(true, [], data);
    var fclTable = $('#flowcltable').DataTable({
      columns: [
        {"data":"1", className: 'firstcol'},
        {"data":"2"}
      ],
      dom: 't',
      data: fclTableData,
      buttons: [
        'copy', 'pdfHtml5','csvHtml5'
      ]
    });
  });
};
