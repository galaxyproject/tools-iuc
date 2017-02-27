// Copyright (c) 2016 Northrop Grumman.
// All rights reserved.

var updateScatter2D = function(){
  scatterData2D['selectedPopulations'] = [];
  scatterDataMFI['selectedPopulations'] = [];
  $('.pop2D').each(function() {
    if (this.checked) {
      scatterData2D['selectedPopulations'].push(parseInt(this.value));
      scatterDataMFI['selectedPopulations'].push(parseInt(this.value));
    }
  });
  processScatterData2D();
  processScatterDataMFI2D();
  displayScatterPlot2D();
};

var processScatterData2D = function() {
  var xData = [],
      yData = [],
      popData = [],
      col1 = [],
      col2 = [],
      pop = [],
      min = Number.MAX_VALUE,
      max = Number.MIN_VALUE;

  min = d3.min(scatterData2D['data'], function(array) {
    return d3.min(array);
  });
  max = d3.max(scatterData2D['data'], function(array) {
    return d3.max(array);
  });
  scatterData2D['min'] = 0;
  scatterData2D['max'] = max;

  col1 = scatterData2D['data'].map(function(value,index) {
    return value[scatterData2D['m1']];
  });
  col2 = scatterData2D['data'].map(function(value,index) {
    return value[scatterData2D['m2']];
  });
  pop = scatterData2D['data'].map(function(value,index) {
    return value[scatterData2D['popCol']];
  });

  for (var i = 0, j = col1.length; i < j; i++) {
    if (scatterData2D['selectedPopulations'].indexOf(pop[i]) >= 0) {
      xData.push(col1[i]);
      yData.push(col2[i]);
      popData.push(pop[i]);
    }
  }

  scatterData2D['popColors'] = popData.map(function(value,index) {
      return color_palette[0][value][0];
  });
  scatterData2D['xData'] = xData;
  scatterData2D['yData'] = yData;
  scatterData2D['popData'] = popData;
  return scatterData2D;
};

var displayScatterToolbar2D = function() {
  $("#xAxisMarker2D").select2();
  $("#yAxisMarker2D").select2();
  $("#view2D").select2();

  scatterData2D['columnHeadings'].map(function(value,index) {
    $('#xAxisMarker2D')
        .append($("<option></option>")
        .attr("value",index)
        .text(value));

    $('#yAxisMarker2D')
        .append($("<option></option>")
        .attr("value",index)
        .text(value));
  });

  $('#xAxisMarker2D').select2("val",0);
  $('#yAxisMarker2D').select2("val",1);

  $("#xAxisMarker2D").on("change",function(e) {
    var m1 = $("#xAxisMarker2D").select2("val");
    scatterData2D['m1'] = m1;
    scatterDataMFI['m1'] = m1;
    updateScatter2D();
  });
  $("#yAxisMarker2D").on("change",function(e) {
    var m2 = $("#yAxisMarker2D").select2("val");
    scatterData2D['m2'] = m2;
    scatterDataMFI['m2'] = m2;
    updateScatter2D();
  });
  $("#view2D").on("change",function(e) {
    var view = $("#view2D").select2("val");
    scatterData2D['view'] = view;
    updateScatter2D();
  });

  $("#updateDisplay2D").on("click",function() {
    $(".pop2D").prop("checked", true);
    $("#selectall2D").prop('checked', true);
    updateScatter2D();
  });
};

var displayScatterPopulation2D = function() {
  $('#populationTable2D tbody').empty();
  scatterData2D['populations'].map(function(value,index) {
    $('#populationTable2D tbody')
        .append('<tr><td align="center">'
            + '<input type="checkbox" checked class="pop2D" value='
            + value + '/></td><td title="' + newNames[value]
            + '">' + newNames[value] + '</td>'
            + '<td><span style="background-color:'
            + color_palette[0][value][0] + '">&nbsp;&nbsp;&nbsp</span></td>'
            + '<td>' + scatterData2D['percent'][value - 1] + '</td></tr>');
  });

  $('#selectall2D').click(function() {
    var checkAll = $("#selectall2D").prop('checked');
    if (checkAll) {
      $(".pop2D").prop("checked", true);
    } else {
      $(".pop2D").prop("checked", false);
    }
    updateScatter2D();
  });

  $('.pop2D').click(function() {
    if ($('.pop2D').length == $(".pop2D:checked").length) {
      $('#selectall2D').prop("checked",true);
    } else {
      $('#selectall2D').prop("checked",false);
    }
    updateScatter2D();
  });

  $('.pop2D').each(function() {
    var selectedpop2D = parseInt(this.value);
    if ($.inArray(selectedpop2D,scatterData2D['selectedPopulations']) > -1) {
      this.checked = true;
    } else {
      this.checked = false;
    }
  });
};

var displayScatterPlot2D = function() {
  var h = $(window).height() - 200,
      w = $("#scatterPlotDiv2D").width(),
      xtitle = scatterData2D['columnHeadings'][scatterData2D['m1']],
      ytitle = scatterData2D['columnHeadings'][scatterData2D['m2']],
      view = scatterData2D['view'],
      layout = {},
      traces = [];

  $("#scatterPlotDiv2D").empty();
  $("#scatterPlotDiv2D").height(h);
  if ( view == 1 || view == 2) {
    var trace1 = {
        x: scatterData2D['xData'],
        y: scatterData2D['yData'],
        mode: 'markers',
        opacity: .75,
        hoverinfo: 'none',
        marker: {
          size: 2,
          color: scatterData2D['popColors']
        },
        type: 'scatter'
    };
    traces.push(trace1);
  }
  if ( view == 1 || view == 3) {
    var trace2 = {
        x: scatterDataMFI['xData'],
        y: scatterDataMFI['yData'],
        mode: 'markers',
        opacity: 1.0,
        hoverinfo: "x+y",
        marker: {
          symbol: 128,
          size: 8,
          color: scatterDataMFI['popColors']
        },
        type: 'scatter'
    };
    traces.push(trace2);
  }

  layout = {
      title: '',
      showlegend: false,
      xaxis: {
          range: [0,scatterData2D['max']],
          title: xtitle
      },
      yaxis: {
          range: [0,scatterData2D['max']],
          title: ytitle
      }
  };

  Plotly.newPlot('scatterPlotDiv2D', traces, layout);
};
