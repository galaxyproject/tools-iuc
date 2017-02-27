// Copyright (c) 2016 Northrop Grumman.
// All rights reserved.

var cl_table = './CLprofiles.txt';
var scores_table = './scores.txt';

var displayCLTable = function(){
  d3.text(cl_table, function(error, data){
    var clHdgs = [],
        clHTML = '<table id="cltable" class="display compact" cellspacing="0" width="100%"/>',
        clTableData = [],
        clHeadings = [];

    if (error){
      alert("Problem retrieving data");
      return;
    }
    clHdgs = data.split("\n")[0].split("\t");
    data = d3.tsv.parse(data);
    clTableData = $.extend(true, [], data);

    clHdgs.forEach(function(d,i){
      clHeadings.push({"data" : d, "title" : d});
    });

    $('#clprofiles').html(clHTML);
    var clTable = $('#cltable').DataTable({
      columns: clHeadings,
      dom: '<"top"Bi>t<"bottom"lp><"clear">',
      pageLength: 25,
      order: [[ 0, "asc" ]],
      data: clTableData,
      buttons: [
        'copy', 'pdfHtml5','csvHtml5'
      ],
      columnDefs: [
        {
          targets: [0,2,3],
          className: "smallcols"
        },
        {
          targets: 4,
          className: "dt-body-left"
        },
        {
          targets: [5,6],
          className: "firstcol"
      }]
    });
  });
};

var displayScoresTable = function(){
  d3.text(scores_table, function(error, data){
    var scoreHTML = '<table id="scoretable" class="display compact" cellspacing="0" width="100%"/>',
        scoreHdgs = [],
        scoreTableData = [],
        scoreHeadings = [];
    if (error){
      alert("Problem retrieving data");
      return;
    }
    scoreHdgs = data.split("\n")[0].split("\t");
    data = d3.tsv.parse(data);

    $('#scores').html(scoreHTML);

    scoreTableData = $.extend(true, [], data);

    scoreHdgs.forEach(function(d,i){
      scoreHeadings.push({"data" : d, "title" : d});
    });

    var scoreTable = $('#scoretable').DataTable({
      columns: scoreHeadings,
      pageLength: 25,
      order: [[ 0, "asc" ]],
      dom: '<"top"Bi>t<"bottom"lp><"clear">',
      data: scoreTableData,
      buttons: [
        'copy', 'pdfHtml5','csvHtml5'
      ],
    });
  });
};
