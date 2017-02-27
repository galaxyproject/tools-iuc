// Copyright (c) 2016 Northrop Grumman.
// All rights reserved.
var updateCSplots = function(plotconfig){
  plotconfig.selectedPopulations = [];
  $(plotconfig.popSelectj).each(function() {
    if (this.checked) {
      plotconfig.selectedPopulations.push(parseInt(this.value));
    }
  });
  if (plotconfig.hasOwnProperty("mtable")) {
    // Update selected markers?
    plotconfig.selectedMarkers = [];
    $(plotconfig.mrkrSelectj).each(function() {
      if (this.checked) {
        plotconfig.selectedMarkers.push(parseInt(this.value));
      }
    });
    // update plot
    updateBoxplot(plotconfig);
  } else {
    updatePlot(plotconfig);
  }
};

var displayPopulationLegend = function(plotconfig) {
  $(plotconfig.table).empty();
  plotconfig.allPopulations.map(function(value,index) {
    $(plotconfig.table).append('<tr><td align="center">'
        + '<input type="checkbox" checked class=' + plotconfig.popSelect
        + ' value=' + value + '/></td><td title="' + newPopNames[value] + '">'
        + newPopNames[value] + '</td><td><span style="background-color:'
        + color_palette[0][value][0] + '">&nbsp;&nbsp;&nbsp;</span></td></tr>');
  });

  $(plotconfig.popSelectAll).click(function() {
    var checkAll = $(plotconfig.popSelectAll).prop('checked');
    if (checkAll) {
      $(plotconfig.popSelectj).prop("checked", true);
    } else {
      $(plotconfig.popSelectj).prop("checked", false);
    }
    updateCSplots(plotconfig);
  });

  $(plotconfig.popSelectj).click(function() {
    if ($(plotconfig.popSelectj).length == $(plotconfig.popSelectCheck).length) {
      $(plotconfig.popSelectAll).prop("checked",true);
    } else {
      $(plotconfig.popSelectAll).prop("checked",false);
    }
    updateCSplots(plotconfig);
  });

  $(plotconfig.popSelectj).each(function() {
    var selectedpopn = parseInt(this.value);
    if ($.inArray(selectedpopn,plotconfig.selectedPopulations) > -1) {
      this.checked = true;
    } else {
      this.checked = false;
    }
  });
};

var displayToolbar = function(plotconfig){
  $(plotconfig.displaybutton).on("click",function() {
    $(plotconfig.popSelectj).prop("checked", true);
    $(plotconfig.popSelectAll).prop("checked", true);
    if (plotconfig.hasOwnProperty("mtable")){
      $(plotconfig.displayMFI).prop("checked", false);
      $(plotconfig.displayvalues).prop("checked", false);
      $(plotconfig.mrkrSelectj).prop("checked", true);
      $(plotconfig.mrkrSelectAll).prop("checked",true);
    }
    updateCSplots(plotconfig);
  });

  if (plotconfig.hasOwnProperty("mtable")){
    $(plotconfig.displayMFI).on("click", function(){
      updateCSplots(plotconfig);
    });
    $(plotconfig.displayvalues).on("click", function(){
      updateCSplots(plotconfig);
    });
  }
  $(plotconfig.toggledisplayj).on("click",function() {
    plotconfig.selectedPopulations = [];
    $(plotconfig.popSelectj).each(function() {
      if (this.checked) {
        plotconfig.selectedPopulations.push(parseInt(this.value));
      }
    })
    if (plotconfig.hasOwnProperty("mtable")){
      plotconfig.selectedMarkers = [];
      $(plotconfig.mrkrSelectj).each(function() {
        if (this.checked) {
          plotconfig.selectedMarkers.push(parseInt(this.value));
        }
      });
      var text = document.getElementById(plotconfig.toggledisplay).firstChild;
      text.data = text.data == "View per marker" ? "View per population" : "View per marker";
      plotconfig.view = plotconfig.view == "p" ? "m" : "p";
      updateBoxplot(plotconfig);
    } else {
      var imgSrc = document.getElementById(plotconfig.toggledisplay);
      imgSrc.src = imgSrc.src.endsWith("stackedsm.png") ? "/static/images/flowtools/barssm.png" : "/static/images/flowtools/stackedsm.png";
      plotconfig.type = plotconfig.type == "barplot" ? "areaplot" : "barplot";
      updatePlot(plotconfig);
    }
  });
  displayPlot(plotconfig);
};

var displayPlot = function(plotconfig) {
  var h = $(window).height() - 200;
  $(plotconfig.plotdivj).empty();
  $(plotconfig.plotdivj).height(h);

  if (plotconfig.hasOwnProperty("mtable")) {
    var nbPop = Object.keys(plotconfig.csdata.mfi[plotconfig.mrkrNames[0]]).length + 2;
    // Get Markers too
    for (var i = 0, nbMarkers = plotconfig.mrkrNames.length; i < nbMarkers; i++) {
      plotconfig.allMarkers.push(i);
      plotconfig.selectedMarkers.push(i);
    }
  } else {
    var nbPop = plotconfig.csdata[0].length;
  }

  for (var i = 2; i < nbPop; i++) {
    plotconfig.allPopulations.push(i - 1);
    plotconfig.selectedPopulations.push(i - 1);
  }

  $(window).on('resize',function() {
    waitForFinalEvent(function() {
      if (plotconfig.hasOwnProperty("mtable")){
          updateBoxplot(plotconfig);
      } else {
          updatePlot(plotconfig);
      }
    }, 500, "resizePlot");
  });

  displayPopulationLegend(plotconfig);
  if (plotconfig.hasOwnProperty("mtable")){
    displayMarkerTable(plotconfig);
    updateBoxplot(plotconfig);
  } else {
    updatePlot(plotconfig);
  }
};

var updatePlot = function(plotconfig) {
  var h = $(window).height() - 200,
      traces = [],
      tmptraces = [],
      x_values = [],
      totals = [];
      layout = {};

  $(plotconfig.plotdivj).empty();
  $(plotconfig.plotdivj).height(h);
  for (var i = 1, j = plotconfig.csdata.length; i < j; i++) {
    x_values.push(newSmpNames[plotconfig.csdata[i][1]]);
  }

  for (var k = 1, i = plotconfig.csdata.length; k < i; k++){
    totals[k] = 0;
    for (var m = 2, o = plotconfig.csdata[0].length; m < o; m++){
      for (var n = 0, p = plotconfig.selectedPopulations.length; n < p; n++){
        if (plotconfig.csdata[0][m] === plotconfig.selectedPopulations[n]) {
          totals[k] += plotconfig.csdata[k][m];
        }
      }
    }
  }

  for (var i = 0, ii = plotconfig.selectedPopulations.length; i < ii; i++) {
    pop = plotconfig.selectedPopulations[i];
    var popName = "Pop " + pop;
    var y_values = [];
    var obj;

    for (var j = 1, jj = plotconfig.csdata.length; j < jj; j++) {
      var newvalue = (plotconfig.csdata[j][pop + 1] / totals[j]) * 100;
      y_values.push(newvalue);
    }
    if (plotconfig.type === "areaplot") {
      obj = {
          x: x_values,
          y: y_values,
          hoverinfo: "x",
          name: popName,
          type: 'area',
          fill: 'tonexty',
          marker: {color: color_palette[0][pop][0]}
      };
    }
    if (plotconfig.type === "barplot") {
      obj = {
          x: x_values,
          y: y_values,
          hoverinfo: "x",
          name: popName,
          type: 'bar',
          marker: {color: color_palette[0][pop][0]}
      };
    }
    tmptraces.push(obj)
  }

  if (plotconfig.type === "barplot") {
    layout = {
        hovermode:'closest',
        title: '',
        barmode: 'stack',
        showlegend: false,
        yaxis: {
            mirror: 'all',
            tickmode: 'array',
            ticktext: ["","20%", "40%", "60%", "80%", "100%"],
            tickvals: [0,20,40,60,80,100],
            title: 'Populations proportions in selected set',
            titlefont: {
                size: 16,
                color: 'grey'
            }
        }
    };
    traces = tmptraces;
  }
  if (plotconfig.type === "areaplot") {
    function stacked(trcs) {
      for(var i=1; i<trcs.length; i++) {
        for(var j=0; j<(Math.min(trcs[i]['y'].length, trcs[i-1]['y'].length)); j++) {
          trcs[i]['y'][j] += trcs[i-1]['y'][j];
        }
      }
      return trcs;
    }
    layout = {
        title: '',
        showlegend: false,
        yaxis: {
            mirror: 'all',
            tickmode: 'array',
            ticktext: ["","20%", "40%", "60%", "80%", "100%"],
            tickvals: [0,20,40,60,80,100],
            title: 'Populations proportions in selected set',
            titlefont: {
                size: 16,
                color: 'grey'
            }
        },
        xaxis: {
            autorange: false,
            range: [-0.2, x_values.length - 0.8]
        }
    };
    traces = stacked(tmptraces);
  }
  Plotly.newPlot(plotconfig.plotdiv,traces,layout);
};
