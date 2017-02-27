// Copyright (c) 2016 Northrop Grumman.
// All rights reserved.
/*
 * Initialize variables for parallelCoordinates display
*/
var pCoordApp = pCoordApp || {};
pCoordApp.allPopulations = [];
pCoordApp.selectedPopulations = [];
pCoordApp.origData;
pCoordApp.flowData;
pCoordApp.headers = [];
pCoordApp.foreground;
pCoordApp.background;
pCoordApp.populations = [];
pCoordApp.allLines;
pCoordApp.lines = [];
pCoordApp.selectedLines = [];

var displayAll = function() {
  displayParallelPlot();
}
/*
 * Display the Population Legend
*/
var displayPopTable = function() {
  $('#popTable tbody').empty();
  pCoordApp.origData.map(function(d,index) {
    $('#popTable tbody')
        .append('<tr><td align="center"><input type="checkbox" '
            + 'id="pop' + d.Population + '" '
            + 'checked class="popSelect" value=' + index + '/></td>'
            + '<td title="' + newNames[d.Population] + '">'
            + newNames[d.Population]
            + '</td><td><span style="background-color:'
            + color_palette[0][index + 1][0]
            + '">&nbsp;&nbsp;&nbsp;</span></td>'
            + '<td>' + d.Percentage + '</td></tr>');
  });

  $('#popSelectAll').click(function() {
    var checkAll = $("#popSelectAll").prop('checked');
    if (checkAll) {
      $(".popSelect").prop("checked", true);
      for (var i = 0; i < pCoordApp.allLines; i ++){
        pCoordApp.selectedLines.push(i);
        pCoordApp.lines.push(i);
      }
    } else {
      $(".popSelect").prop("checked", false);
      pCoordApp.selectedLines = [];
      pCoordApp.lines = [];
    }

    pCoordApp.selectedPopulations = [];
    $('.popSelect').each(function() {
      if (this.checked) {
        pCoordApp.selectedPopulations.push(parseInt(this.value));
      }
    });
    displayTableGrid();
    if (checkAll) {
      displayParallelPlot();
    } else {
      updateParallelForeground();
    }
  });

  $('.popSelect').click(function() {
    if ($('.popSelect').length == $(".popSelect:checked").length) {
        $('#popSelectAll').prop("checked",true);
    } else {
        $('#popSelectAll').prop("checked",false);
    }

    pCoordApp.selectedPopulations = [];
    $('.popSelect').each(function() {
      if (this.checked) {
        pCoordApp.selectedPopulations.push(parseInt(this.value));
      }
    });

    pCoordApp.selectedLines = [];
    pCoordApp.lines = [];

    pCoordApp.origData.forEach(function(d,idx){
      if ($.inArray(pCoordApp.populations.indexOf(d.Population), pCoordApp.selectedPopulations) > -1) {
        pCoordApp.selectedLines.push(idx);
        pCoordApp.lines.push(idx);
      }
    });

    displayTableGrid();
    updateParallelForeground();
  });
  updatePopTable();
};

var updatePopTable = function() {
  $('.popSelect').each(function() {
    var pop = parseInt(this.value),
        selectedPops = pCoordApp.origData.map(function(d){
      if ($.inArray(d.idx, pCoordApp.selectedLines) > -1){
        return pCoordApp.populations.indexOf(d.Population);
      }
    });
    if ($.inArray(pop,selectedPops) > -1) {
      this.checked = true;
    } else {
      this.checked = false;
    }
  });
};

/*
 * Display the table under the graph
*/
var displayTableGrid = function() {
  var updatedData = [],
      displayData = [],
      colNames = [],
      pctargets = [],
      colTable = [],
      tableHTML = [],
      textCol = [],
      colOrder = [],
      targetCol = 0;

  $("#tableDiv").empty();
  updatedData = $.extend(true, [], tableContent);
  updatedData.forEach(function(d, idx){d.idx = idx});
  displayData = updatedData.filter(function(d, index) {
    if ($.inArray(index,pCoordApp.selectedLines) > -1) {
      return d;
    }
  });

  targetCol = pCoordApp.headers.length - 2;
  pCoordApp.headers.forEach(function(d,i){
    colTable.push("<th>" + d + "</th>");
    colNames.push({"data":d});
    if (i < targetCol){
      pctargets.push(i);
    }
  });
  textCol = [targetCol, targetCol + 1];
  colOrder = textCol.concat(pctargets);
  tableHTML = [
      '<table id="pcTable" class="pctable display compact" cellspacing="0" width="100%">',
      '<thead>',
      '<tr>',
      colTable.join("\n"),
      '</tr>',
      '</thead>',
      '</table>',
  ];

  $('#tableDiv').html(tableHTML.join("\n"));
  var pcTable = $('#pcTable').DataTable({
      columns: colNames,
      data: displayData,
      order: [[ targetCol, "asc" ]],
      pageLength: 10,
      //paging: false,
      scrollY: 250,
      scrollCollapse: true,
      scrollX: true,
      dom: '<"top"B>t<"bottom"lip><"clear">',
      columnDefs: [{
          targets: pctargets,
          className: "dt-body-right",
        }, {
          targets: [targetCol, targetCol+1],
          className: "dt-body-center"
      }],
      buttons: [
          'copy', 'pdfHtml5','csvHtml5', 'colvis'
      ],
      colReorder: {
          order:colOrder
      },
      select: true
  });

  $('#pcTable').on('mouseover', 'tr', function() {
    var data = pcTable.row(this).data();
    if (data != undefined) {
      var line = data.idx;
      pCoordApp.selectedLines = [ line ];
      updateParallelForeground();
    }
  });
  $('#pcTable').on('mouseleave', 'tr', function() {
    pCoordApp.selectedLines = [];
    for (var i = 0, j = pCoordApp.lines.length; i < j; i++) {
      pCoordApp.selectedLines.push(pCoordApp.lines[i]);
    }
    updateParallelForeground();
  });
};

/*
 * Display The Main Plot
*/
var displayParallelPlot = function() {
  var margin = {top: 30, right: 10, bottom: 10, left: 10},
      h = $("#chartDiv").height()/1.5,
      w = $("#plotDiv").width(),
      width = w - margin.left - margin.right,
      height = h - margin.top - margin.bottom,
      dragging = {},
      y = {};

  $("#plotDiv").empty();
  $("#plotDiv").height(h);
  var svg = d3.select("#plotDiv").append("svg")
      .attr("width", width + margin.left + margin.right)
      .attr("height", height + margin.top + margin.bottom)
    .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  // Y axis label
  svg.append("text")
      .attr("class", "ylabel")
      .attr("transform", "rotate(-90)")
      .attr("y", 0 - margin.left)
      .attr("x", 0 - (height / 2))
      .attr("dy", "1em")
      .style("text-anchor", "middle")
      .text("MFI");

  var x = d3.scale.ordinal().rangePoints([0, width], 1);

  // Use this to scale line width to percentage population
  var pd = d3.extent(pCoordApp.origData, function(p) {
    return +p['Percentage'];
  });
  var popScale = d3.scale.linear().range([1,5]).domain(pd);

  var line = d3.svg.line();
  var axis = d3.svg.axis().orient("left").ticks(8);

  var dimensions = d3.keys(pCoordApp.flowData[0]).filter(function(d) {
    return (y[d] = d3.scale.linear()
        .domain(d3.extent(pCoordApp.flowData,function(p) { return +p[d]; }))
        .range([height, 0]));
  });
  x.domain(dimensions);

  function path(d) {
    return line(dimensions.map(function(p) {
      return [x(p), y[p](d[p])];
    }));
  }
  function position(d) {
    var v = dragging[d];
    return v == null ? x(d) : v;
  }
  function transition(g) {
    return g.transition().duration(500);
  }
  function brush() {
    var actives = dimensions.filter(function(p) {
      return !y[p].brush.empty();
    });
    var extents = actives.map(function(p) {
      return y[p].brush.extent();
    });
    var indices = pCoordApp.origData.filter(function(d) {
      var line = d.idx;
      var tf = actives.every(function(p,i) {
        return extents[i][0] <= pCoordApp.flowData[line][p] &&
                          pCoordApp.flowData[line][p] <= extents[i][1];
      });
      if (tf) {
        return line.toString();
      }
    });
    pCoordApp.selectedLines = indices.map(function(d) {
      return d.idx;
    });
    pCoordApp.lines = indices.map(function(d) {
      return d.idx;
    });

    updateParallelForeground();
    updatePopTable();
    displayTableGrid();
  };

  // Display paths in light gray color, to use as reference
  pCoordApp.background = svg.append("g")
      .attr("class", "background")
    .selectAll("path")
      .data(pCoordApp.flowData)
    .enter().append("path")
      .attr("d", path);

  // Add foreground lines for focus, color by population.
  pCoordApp.foreground = svg.append("g")
      .attr("class", "foreground")
    .selectAll("path")
      .data(pCoordApp.origData)
    .enter().append("path")
      .attr("d", path)
      .attr("stroke",function(d){
        var pop = pCoordApp.populations.indexOf(d.Population) + 1;
        return color_palette[0][pop][0];
      })
      //.attr("stroke-width", 2);
      // Use this if you want to scale the lines based on
      // population percentage
      .attr("stroke-width", function(d) {
        var pop = pCoordApp.populations.indexOf(d.Population);
        var w = popScale(pCoordApp.origData[pop]['Percentage']);
        w = parseInt(w);
        return w;
      });

  // Add a group element for each dimension.
  var g = svg.selectAll(".dimension")
      .data(dimensions)
    .enter().append("g")
      .attr("class", "dimension")
      .attr("transform", function(d) { return "translate(" + x(d) + ")"; })
      .call(d3.behavior.drag()
      .origin(function(d) { return {x: x(d)}; })
      .on("dragstart", function(d) {
        dragging[d] = x(d);
        pCoordApp.background.attr("visibility", "hidden");
      })
      .on("drag", function(d) {
        dragging[d] = Math.min(width, Math.max(0, d3.event.x));
        pCoordApp.foreground.attr("d", path);
        dimensions.sort(function(a, b) {
          return position(a) - position(b);
        });
        x.domain(dimensions);
        g.attr("transform", function(d) {
          return "translate(" + position(d) + ")";
        });
      })
      .on("dragend", function(d) {
        delete dragging[d];
        transition(d3.select(this))
            .attr("transform", "translate(" + x(d) + ")");
        transition(pCoordApp.foreground)
            .attr("d", path);
        pCoordApp.background.attr("d", path)
          .transition()
            .delay(500)
            .duration(0)
            .attr("visibility", null);
    }));

  // Add an axis and title.
  g.append("g")
      .attr("class", "axis")
      .each(function(d) {
        d3.select(this).call(axis.scale(y[d]));
      });
  g.append("g")
      .attr("class", "xlabel")
    .append("text")
      .style("text-anchor", "middle")
      .attr("y", -9)
      .text(function(d) { return d; });

  // Add and store a brush for each axis.
  g.append("g")
      .attr("class", "brush")
      .each(function(d) { d3.select(this).call(y[d].brush = d3.svg.brush().y(y[d]).on("brush", brush)); })
    .selectAll("rect")
      .attr("x", -8)
      .attr("width", 16);

  // Control line opacity.
  $('#pcline_opacity').on('change', (function() {
    var val = $(this).val();
    $('#plotDiv .foreground path').css('stroke-opacity', val.toString());
    $('#pcopacity').html((Math.round(val*10000)/100) + "%");
  }));
};

var updateParallelForeground = function() {
  pCoordApp.foreground[0].map(function(d) {
    var ln = parseInt(d['__data__']['idx']);
    if ($.inArray(ln, pCoordApp.selectedLines) < 0){
      d.style.display = "none";
    } else {
      d.style.display = null;
    }
  });
};

/*
 * Retrieve the data, then call display functions
*/
var displayParallelCoordinates = function() {
  pCoordApp.origData = $.extend(true,[], tableContent);
  pCoordApp.headers = Object.keys(pCoordApp.origData[0]);
  pCoordApp.origData.forEach(function(d,idx) {
    d.idx = idx;
    pCoordApp.selectedLines.push(idx);
    pCoordApp.lines.push(idx);
    if (!pCoordApp.populations.includes(d.Population)){
      pCoordApp.populations.push(d.Population);
    }
  });
  /*
   * For the plot use only the MFI information
   * for each populations. Store in flowData
  */
  pCoordApp.flowData = $.extend(true,[],tableContent);
  pCoordApp.flowData.forEach(function(d, idx) {
    delete d['Population'];
    delete d['Count'];
    delete d['Percentage'];
    delete d.Comment;
    pCoordApp.allPopulations.push(idx);
    pCoordApp.selectedPopulations.push(idx);
    pCoordApp.selectedLines.push(idx);
    pCoordApp.lines.push(idx);
  });

  pCoordApp.allLines = pCoordApp.flowData.length;
  displayPopTable();
  displayTableGrid();
  displayParallelPlot();

  $("#resetPCoordDisplay").on("click",function() {
    for (var i = 0; i < pCoordApp.allLines; i++) {
      pCoordApp.allPopulations.push(i);
      pCoordApp.selectedPopulations.push(i);
      pCoordApp.selectedLines.push(i);
      pCoordApp.lines.push(i);
    }
    $("#popSelectAll").prop('checked',true);
    $(".popSelect").prop("checked",true);

    var opcty = ".8";
    $('#plotDiv .foreground path').css('stroke-opacity', opcty);
    $('#pcopacity').html("80%");
    $('#pcline_opacity').val(0.8);

    displayPopTable();
    displayTableGrid();
    displayParallelPlot();
  });

  $(window).on('resize',function() {
    waitForFinalEvent(function() {
      displayAll();
    }, 500, "resizeParallelCoordinates");
  });
}
