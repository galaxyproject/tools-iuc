// Copyright (c) 2016 Northrop Grumman.
// All rights reserved.

var updateScatter3D = function(){
  scatterData3D['selectedPopulations'] = [];
  scatterData3DMFI['selectedPopulations'] = [];
  $('.pop3D').each(function() {
    if (this.checked) {
      scatterData3D['selectedPopulations'].push(parseInt(this.value));
      scatterData3DMFI['selectedPopulations'].push(parseInt(this.value));
    }
  });
  processScatterData3D();
  processScatterData3DMFI();
  displayScatterPlot3D();
};

var processScatterData3D = function() {
  var xData = [],
      yData = [],
      zData = [],
      popData = [],
      pop = [],
      col1 = [],
      col2 = [],
      col3 = [],
      min = Number.MAX_VALUE,
      max = Number.MIN_VALUE;

  min = d3.min(scatterData3D['data'], function(array) {
    return d3.min(array);
  });
  max = d3.max(scatterData3D['data'], function(array) {
    return d3.max(array);
  });
  scatterData3D['min'] = 0;
  scatterData3D['max'] = max;

  col1 = scatterData3D['data'].map(function(value,index) {
    return value[scatterData3D['m1']];
  });
  col2 = scatterData3D['data'].map(function(value,index) {
    return value[scatterData3D['m2']];
  });
  col3 = scatterData3D['data'].map(function(value,index) {
    return value[scatterData3D['m3']];
  });
  pop = scatterData3D['data'].map(function(value,index) {
    return value[scatterData3D['popCol']];
  });

  for (var i = 0, j = col1.length; i < j; i++) {
    if (scatterData3D['selectedPopulations'].indexOf(pop[i]) >= 0) {
      xData.push(col1[i]);
      yData.push(col2[i]);
      zData.push(col3[i]);
      popData.push(pop[i]);
    }
  }

  scatterData3D['popColors'] = popData.map(function(value,index) {
    return color_palette[0][value][0];
  });
  scatterData3D['xData'] = xData;
  scatterData3D['yData'] = yData;
  scatterData3D['zData'] = zData;
  scatterData3D['popData'] = popData;
  return scatterData3D;
};

var displayScatterToolbar3D = function() {
  $("#xAxisMarker3D").select2();
  $("#yAxisMarker3D").select2();
  $("#zAxisMarker3D").select2();
  $("#view3D").select2();

  scatterData3D['columnHeadings'].map(function(value,index) {
    $('#xAxisMarker3D')
        .append($("<option></option>")
        .attr("value",index)
        .text(value));

    $('#yAxisMarker3D')
        .append($("<option></option>")
        .attr("value",index)
        .text(value));

    $('#zAxisMarker3D')
        .append($("<option></option>")
        .attr("value",index)
        .text(value));
  });

  $('#xAxisMarker3D').select2("val",0);
  $('#yAxisMarker3D').select2("val",1);
  $('#zAxisMarker3D').select2("val",2);

  $("#xAxisMarker3D").on("change",function(e) {
    var m1 = $("#xAxisMarker3D").select2("val");
    scatterData3D['m1'] = m1;
    scatterData3DMFI['m1'] = m1;
    updateScatter3D();
  });
  $("#yAxisMarker3D").on("change",function(e) {
    var m2 = $("#yAxisMarker3D").select2("val");
    scatterData3D['m2'] = m2;
    scatterData3DMFI['m2'] = m2;
    updateScatter3D();
  });
  $("#zAxisMarker3D").on("change",function(e) {
    var m3 = $("#zAxisMarker3D").select2("val");
    scatterData3D['m3'] = m3;
    scatterData3DMFI['m3'] = m3;
    updateScatter3D();
  });

  $("#view3D").on("change",function(e) {
    var view = $("#view3D").select2("val");
    scatterData3D['view'] = view;
    updateScatter3D();
  });

  $("#updateDisplay3D").on("click",function() {
    $(".pop3D").prop("checked", true);
    $("#selectall3D").prop('checked', true);
    updateScatter3D();
  });
};

var displayScatterPopulation3D = function() {
  $("#populationTable3D tbody").empty();
  scatterData3D['populations'].map(function(value) {
    $('#populationTable3D tbody')
        .append('<tr><td align="center">'
                + '<input type="checkbox" checked class="pop3D" value='
                + value + '/></td><td title="'+ newNames[value]
                + '">'+ newNames[value] + '</td>'
                + '<td><span style="background-color:'
                + color_palette[0][value][0] + '">&nbsp;&nbsp;&nbsp</span></td>'
                + '<td>' + scatterData3D['percent'][value - 1] + '</td></tr>');
  });

  $('#selectall3D').click(function() {
    var checkAll = $("#selectall3D").prop('checked');
    if (checkAll) {
      $(".pop3D").prop("checked", true);
    } else {
      $(".pop3D").prop("checked", false);
    }
    updateScatter3D();
  });
  $('.pop3D').click(function() {
    if ($('.pop3D').length == $(".pop3D:checked").length) {
      $('#selectall3D').prop("checked",true);
    } else {
      $('#selectall3D').prop("checked",false);
    }
    updateScatter3D();
  });

  $('.pop3D').each(function() {
    var selectedpop3D = parseInt(this.value);
    if ($.inArray(selectedpop3D,scatterData3D['selectedPopulations']) > -1) {
      this.checked = true;
    } else {
      this.checked = false;
    }
  });
};

var displayScatterPlot3D = function() {
  var h = $(window).height() - 200,
      w = $("#scatterPlotDiv3D").width(),
      xtitle = scatterData3D['columnHeadings'][scatterData3D['m1']],
      ytitle = scatterData3D['columnHeadings'][scatterData3D['m2']],
      ztitle = scatterData3D['columnHeadings'][scatterData3D['m3']],
      view = scatterData3D['view'],
      traces = [],
      layout = {};

  $("#scatterPlotDiv3D").empty();
  $("#scatterPlotDiv3D").height(h);

  if ( view == 1 || view == 2 ) {
    var trace1 = {
        x: scatterData3D['xData'],
        y: scatterData3D['yData'],
        z: scatterData3D['zData'],
        mode: 'markers',
        opacity: .75,
        hoverinfo: "none",
        marker: {
          size: 2,
          color: scatterData3D['popColors']
        },
        type: 'scatter3d'
    };
    traces.push(trace1);
  };

  if ( view == 1 || view == 3) {
    var trace2 = {
        x: scatterData3DMFI['xData'],
        y: scatterData3DMFI['yData'],
        z: scatterData3DMFI['zData'],
        mode: 'markers',
        opacity: 1.0,
        hoverinfo: "x+y+z",
        marker: {
          symbol: "circle-open",
          size: 8,
          color: scatterData3DMFI['popColors']
        },
        type: 'scatter3d'
    };
    traces.push(trace2);
  }

  layout = {
      title: '',
      showlegend: false,
      scene: {
        aspectratio: {
          x: 1,
          y: 1,
          z: 1
        },
        camera: {
          center: {
            x: 0,
            y: 0,
            z: 0
          },
          eye: {
            x: 1.25,
            y: 1.25,
            z: 1.25
          },
          up: {
            x: 0,
            y: 0,
            z: 1
          }
        },
        xaxis: {
          type: 'linear',
          title: xtitle,
          range: [0, scatterData3D['max']],
          zeroline: false
        },
        yaxis: {
          type: 'linear',
          title: ytitle,
          range: [0, scatterData3D['max']],
          zeroline: false
        },
        zaxis: {
          type: 'linear',
          title: ztitle,
          range: [0, scatterData3D['max']],
          zeroline: false
        }
      }
  };

  Plotly.newPlot('scatterPlotDiv3D', traces, layout);
};
