// Copyright (c) 2016 Northrop Grumman.
// All rights reserved.

var scatterData2D = {};
var scatterData3D = {};
var scatterData3DMFI = {};
var scatterDataMFI = {};
var tableContent;
var newNames ={};
var configbp = {};
var bpurl = "./boxplotData.json";

var waitForFinalEvent = (function () {
  var timers = {};
  return function (callback, ms, uniqueId) {
    if (!uniqueId) {
      uniqueId = "Don't call this twice without a uniqueId";
    }
    if (timers[uniqueId]) {
      clearTimeout (timers[uniqueId]);
    }
    timers[uniqueId] = setTimeout(callback, ms);
  };
})();

var updateOverviewPlotDisplay = function(data) {
  $('#overviewPlotDiv').empty();
  $('#overviewPlotDiv').html(data);
};

var displayMFI = function() {
  var url = "flow.mfi_pop";
  d3.text(url, function(error, data){
    var mfiHdgs = [],
        pp = [],
        mfiTableData = [],
        mfiTableHeadings = [],
        mfiTargets = [],
        mfiTableHTML = '<table id="mfitable" class="dtable display compact" cellspacing="0" width="100%"/>',
        popcol = 1,
        mfiEditorData =[];

    if (error){
      alert("Problem retrieving data");
      return;
    }
    mfiHdgs = data.split("\n")[0].split("\t");
    pp = mfiHdgs.pop();
    mfiHdgs.unshift(pp);
    mfiHdgs.unshift("Comment");
    data = d3.tsv.parse(data);
    function handleSubmit(method, url, d, successCallBack, errorCallBack) {
      var output = {data : mfiTableData};
      successCallBack(output);
    }

    $('#mfiDiv').empty();
    mfiTableData = $.extend(true,[],data);
    mfiTableData.forEach(function(d){
      d.Comment = d.Population;
      newNames[parseInt(d.Population)] = d.Comment;
    });
    //var mfiData = mfiTableData.filter(function(d){return d});
    tableContent = $.extend(true, [], mfiTableData);
    mfiHdgs.forEach(function(d,i) {
      mfiTableHeadings.push({"data": d, "title": d});
      mfiEditorData.push({"label" : d, "name" : d});
    });

    for (var i = 2, j = mfiHdgs.length - 2; i<j; i++){
      mfiTargets.push(i);
    }

    $('#mfiDiv').html(mfiTableHTML);
    var editor = new $.fn.dataTable.Editor({
        ajax: handleSubmit,
        table: '#mfitable',
        fields: mfiEditorData,
        idSrc: 'Population'
    });
    $('#mfitable').on('click', 'tbody td:first-child', function (e) {
      editor.bubble(this);
    });
    var mfiTable = $('#mfitable').DataTable({
        columns: mfiTableHeadings,
        data: mfiTableData,
        order: [[ popcol, "asc" ]],
        pageLength: 25,
        dom: '<"top"Bi>t<"bottom"lp><"clear">',
        columnDefs: [{
            targets: mfiTargets,
            className: "dt-body-right",
            render: function(data, type, row){
              return parseFloat(data).toFixed(2);
            }
          }, {
            targets: [mfiHdgs.length - 2],
            className: "dt-body-right"
          }, {
            targets: [mfiHdgs.length-1],
            className: "dt-body-right",
            render: function(data, type, row){
                    return parseFloat(data).toFixed(2) + '%';
            }
          }
        ],
        buttons: [
            'copy', 'pdfHtml5','csvHtml5', 'colvis'
        ],
        colReorder: {fixedColumnsLeft:1},
        select: true
    });
    editor.on('preSubmit', function(e, object, action){
      var data = object.data;
      var key = Object.keys(data)[0];
      var count = object.data[key]['Comment'];
      mfiTableData.forEach(function(d) {
        if (d.Population === key) {
          d.Comment = count;
          newNames[parseInt(d.Population)] = count;
        }
      });
      tableContent = $.extend(true, [], mfiTableData);
    });
  });
};

var displayOverviewPlot = function() {
  var url = "flow.overview";
  $.ajax({
    url: url,
    dataType: "text",
    success: function(data) {
      updateOverviewPlotDisplay(data);
    }
  });
};

var displayScatter2D = function() {
  var url = "flow.sample";
  $.ajax({
    url: url,
    dataType: "text",
    success: function(text) {
      var mfi_url = "flow.mfi_pop";
      $.ajax({
        url: mfi_url,
        dataType: "text",
        success: function(mfi_text) {
          scatterDataMFI = new processMFI(mfi_text);
          scatterData2D = new processData(text);
          displayScatterToolbar2D();
          displayScatterPopulation2D();
          processScatterData2D();
          processScatterDataMFI2D();
          displayScatterPlot2D();
          $(window).on('resize',function() {
            waitForFinalEvent(function() {
              processScatterData2D();
              displayScatterPlot2D();
            }, 500, "resize2D");
          });
        }
      });
    }
  });
};

var displayScatter3D = function() {
  var url = "flow.sample";
  $.ajax({
    url: url,
    dataType: "text",
    success: function(text) {
      var mfi_url = "flow.mfi_pop";
      $.ajax({
        url: mfi_url,
        dataType: "text",
        success: function(mfi_text) {
          scatterData3DMFI = new processMFI(mfi_text);
          scatterData3D = new processData(text);
          displayScatterToolbar3D();
          displayScatterPopulation3D();
          processScatterData3D();
          processScatterData3DMFI();
          displayScatterPlot3D();
          $(window).on('resize',function() {
            waitForFinalEvent(function() {
              processScatterData3D();
              displayScatterPlot3D();
            }, 500, "resize3D");
          });
        }
      });
    }
  });
};

var displayMFIBoxplot = function() {
  $.ajax({
    url: bpurl,
    dataType: "json",
    success: function(data) {
      configbp = {
        displaybutton : '#updateDisplaybp',
        toggledisplayj : '#changeDisplay',
        toggledisplay : 'changeDisplay',
        popSelectj : '.popSelectbp',
        plotdivj : '#plotDivbp',
        csdata : data,
        plotdiv : 'plotDivbp',
        type : 'boxplot',
        table : '#popTablebp tbody',
        popSelect : 'popSelectbp',
        allMarkers : [],
        selectedMarkers: [],
        allPopulations : [],
        selectedPopulations : [],
        popSelectAll : '#popSelectAllbp',
        popSelectCheck: '.popSelectbp:checked',
        mrkrSelectAll : '#mrkrSelectAll',
        mrkrSelectCheck: '.mrkrSelect:checked',
        mrkrSelect : 'mrkrSelect',
        mtable : '#mrkrTable tbody',
        mrkrSelectj: '.mrkrSelect',
        displayvalues: '#displayLabels',
        displayMFI: '#displayMFI',
        view: 'p',
        mrkrNames :  Object.keys(data.mfi)
      };
      displayToolbar(configbp);
    }
  });
};

function processData(text) {
  var data = d3.tsv.parseRows(text).map(function(row) {
    return row.map(function(value) {
      if (isNaN(value)) {
        return value;
      }
      return +value;
    });
  });

  this.columnHeadings = data.shift();
  this.columnHeadings.pop();
  var popCol = data[0].length - 1;
  var p = data.map(function(value,index) {
    return parseInt(value[popCol]);
  });

  var populations = {};
  for (var i = 0; i < p.length; i++) {
    if (populations[p[i]] === undefined) {
      populations[p[i]] = 1;
    } else {
      populations[p[i]] = populations[p[i]] + 1;
    }
  }

  this.popCol = popCol;
  this.populations = d3.set(p).values();
  this.populations = this.populations.map(function(value,index) {
    return parseInt(value);
  });
  this.selectedPopulations = this.populations;
  this.percent = this.populations.map(function(value,index) {
    return Math.floor(populations[value] * 10000.0 / data.length) / 100.0;
  });

  this.data = data;
  this.m1 = 0;
  this.m2 = 1;
  this.m3 = 2;
  this.view = 1;
};

function processMFI(text) {
  data = d3.tsv.parseRows(text).map(function(row) {
    return row.map(function(value) {
      if (isNaN(value)) {
        return value;
      }
      return +value;
    });
  });

  // Get the Headings Row, then remove the Count, Percentage and
  // Population headings
  this.columnHeadings = data.shift();
  this.columnHeadings.pop();
  this.columnHeadings.pop();
  this.columnHeadings.pop();

  var popCol = data[0].length -1;
  var pop = data.map(function(value,index) {
    return parseInt(value[popCol]);
  });

  var perCol = data[0].length -2;
  var per = data.map(function(value,index) {
    return parseFloat(value[perCol]);
  });

  var countCol = data[0].length -3;
  var count = data.map(function(value,index) {
    return parseInt(value[countCol]);
  });

  this.popCol = popCol;
  this.populations = pop;
  this.selectedPopulations = pop;
  this.percent = per;
  this.counts = count;

  var l = data[0].length;
  this.data = data.map(function(row) {
    return row.splice(0,countCol);
  });
  this.poplist = pop;
  this.m1 = 0;
  this.m2 = 1;
  this.m3 = 2;
};
