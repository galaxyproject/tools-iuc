// Copyright (c) 2016 Northrop Grumman.
// All rights reserved.
/*
 * Initialize variables for parallelCoordinates display
*/
var pcAppMFI = pcAppMFI || {};

pcAppMFI.origData;
pcAppMFI.flowData;
pcAppMFI.background;
pcAppMFI.foreground;
pcAppMFI.selectedLines = [];
pcAppMFI.selectedPopulations = [];
pcAppMFI.selectedSamples = [];
pcAppMFI.populations = [];
pcAppMFI.samples = [];
pcAppMFI.lines = [];
pcAppMFI.allLines;
pcAppMFI.headers = [];

var displayAllm = function() {
  displayParallelPlotm();
};
/*
 * Display the Population Legend
*/
var displayPopTablem = function() {
  $('#popTablePCm tbody').empty();
  pcAppMFI.populations.map(function(d, index) {
    $('#popTablePCm tbody').append('<tr><td align="center">'
        + '<input type="checkbox" '
        + 'id="'+ d + '" '
        + 'checked class="popSelectPCm" value='
        + index + '/></td><td title="' + newPopNames[d]
        + '">' + newPopNames[d]
        + '</td><td><span style="background-color:'
        + color_palette[0][index + 1][0]
        + '">&nbsp;&nbsp;&nbsp;</span></td></tr>');
  });

  $('#popSelectAllPCm').click(function() {
    var checkAll = $("#popSelectAllPCm").prop('checked');
    if (checkAll) {
      $(".popSelectPCm").prop("checked", true);
      for (var i = 0; i < pcAppMFI.allLines; i ++) {
        pcAppMFI.selectedLines.push(i);
        pcAppMFI.lines.push(i);
      }
    } else {
      $(".popSelectPCm").prop("checked", false);
      pcAppMFI.selectedLines = [];
      pcAppMFI.lines = [];
    }

    pcAppMFI.selectedPopulations = [];
    $('.popSelectPCm').each(function() {
      if (this.checked) {
        pcAppMFI.selectedPopulations.push(parseInt(this.value));
      }
    });

    displayTableGridm();
    if (checkAll) {
      displayParallelPlotm();
    } else {
      updateParallelForegroundidx();
    }
  });

  $('.popSelectPCm').click(function() {
    if ($('.popSelectPCm').length == $(".popSelectPCm:checked").length) {
      $('#popSelectAllPCm').prop("checked",true);
    } else {
      $('#popSelectAllPCm').prop("checked",false);
    }
    pcAppMFI.selectedPopulations = [];
    $('.popSelectPCm').each(function() {
      if (this.checked) {
        pcAppMFI.selectedPopulations.push(parseInt(this.value));
      }
    });
    pcAppMFI.selectedLines = [];
    pcAppMFI.lines = [];
    pcAppMFI.origData.forEach(function(d,idx){
      if ($.inArray(pcAppMFI.populations.indexOf(d.Population), pcAppMFI.selectedPopulations) > -1) {
        if ($.inArray(pcAppMFI.samples.indexOf(d.SmpName), pcAppMFI.selectedSamples) > -1){
          pcAppMFI.selectedLines.push(idx);
          pcAppMFI.lines.push(idx);
        }
      }
    });
    displayTableGridm();
    updateParallelForegroundidx();
  });
  updatePopTableidx();
  updateSmpTableidx();
};

var updatePopTableidx = function() {
  $('.popSelectPCm').each(function() {
    var pop = parseInt(this.value);
    var selectedPops = pcAppMFI.origData.map(function(d){
      if ($.inArray(d.idx, pcAppMFI.selectedLines) > -1){
        return pcAppMFI.populations.indexOf(d.Population);
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
* Display Sample Legend
*/
var displaySmpTablem = function(){
  $('#smpTablePCm tbody').empty();
  pcAppMFI.samples.map(function(d, index) {
    $('#smpTablePCm tbody').append('<tr><td title="'
        + newSmpNames[d] + '">' + newSmpNames[d]
        + '</td><td align="center">' + '<input type="checkbox" '
        + 'id="' + d + '" ' + 'checked class="smpSelectPCm" value='
        + index + '></td></tr>');
  });

  $('#smpSelectAllPCm').click(function() {
    var checkAll = $("#smpSelectAllPCm").prop('checked');
    if (checkAll) {
      $(".smpSelectPCm").prop("checked", true);
      for (var i = 0; i < pcAppMFI.allLines; i ++) {
        pcAppMFI.selectedLines.push(i);
        pcAppMFI.lines.push(i);
      }
    } else {
      $(".smpSelectPCm").prop("checked", false);
      pcAppMFI.selectedLines = [];
      pcAppMFI.lines = [];
    }
    pcAppMFI.selectedSamples = [];
    $('.smpSelectPCm').each(function() {
      if (this.checked) {
        pcAppMFI.selectedSamples.push(parseInt(this.value));
      }
    });
    displayTableGridm();
    if (checkAll) {
      displayParallelPlotm();
    } else {
      updateParallelForegroundidx();
    }
  });

  $('.smpSelectPCm').click(function() {
    if ($('.smpSelectPCm').length == $(".smpSelectPCm:checked").length) {
      $('#smpSelectAllPCm').prop("checked",true);
    } else {
      $('#smpSelectAllPCm').prop("checked",false);
    }
    pcAppMFI.selectedSamples = [];
    $('.smpSelectPCm').each(function() {
      if (this.checked) {
        pcAppMFI.selectedSamples.push(parseInt(this.value));
      }
    });
    pcAppMFI.selectedLines = [];
    pcAppMFI.lines = [];
    pcAppMFI.origData.forEach(function(d,idx) {
      if ($.inArray(pcAppMFI.populations.indexOf(d.Population), pcAppMFI.selectedPopulations) > -1) {
        if ($.inArray(pcAppMFI.samples.indexOf(d.SmpName), pcAppMFI.selectedSamples) > -1){
          pcAppMFI.selectedLines.push(idx);
          pcAppMFI.lines.push(idx);
        }
      }
    });
    displayTableGridm();
    updateParallelForegroundidx();
  });
};

var updateSmpTableidx = function() {
  $('.smpSelectPCm').each(function() {
    var smp = parseInt(this.value),
        selectedSamples = pcAppMFI.origData.map(function(d){
      if ($.inArray(d.idx, pcAppMFI.selectedLines) > -1){
        return pcAppMFI.samples.indexOf(d.SmpName);
      }
    });
    if ($.inArray(smp,selectedSamples) > -1) {
      this.checked = true;
    } else {
      this.checked = false;
    }
  });
};
/*
 * Display Data table
*/
var displayTableGridm = function() {
  var colTablem = [],
      colNamesm = [],
      pctargetsm = [],
      displayDatamfi = [],
      targetColm = pcAppMFI.headers.length - 3,
      textColm = [],
      colOrderm = [],
      tableHTMLm = [];

  $("#tableDivPCm").empty();

  displayDatamfi = pcAppMFI.origData.filter(function(d,i) {
    if ($.inArray(i,pcAppMFI.selectedLines) > -1) {
      return d;
    }
  });
  displayDatamfi.forEach(function(d){
    d.EditedPopName = newPopNames[d.Population];
    d.SampleName = newSmpNames[d.SmpName];
  });
  pcAppMFI.headers.forEach(function(d,i){
    colTablem.push("<th>" + d + "</th>");
    colNamesm.push({"data":d});
    if (i < targetColm - 1){
      pctargetsm.push(i);
    }
  });
  textColm = [targetColm, targetColm + 1, targetColm + 2];
  colOrderm = textColm.concat(pctargetsm);
  colOrderm.push(targetColm - 1);
  tableHTMLm = [
      '<table id="pcTableMFI" class="pctable display compact nowrap" cellspacing="0" width="100%">',
      '<thead>',
      '<tr>',
      colTablem.join("\n"),
      '</tr>',
      '</thead>',
      '</table>',
  ];
  $('#tableDivPCm').html(tableHTMLm.join("\n"));
  var pcTablem = $('#pcTableMFI').DataTable({
      columns: colNamesm,
      data: displayDatamfi,
      order: [[ targetColm, "asc" ]],
      pageLength: 10,
      //paging: false,
      scrollY: 250,
      scrollCollapse: true,
      scrollX: true,
      dom: '<"top"B>t<"bottom"lip><"clear">',
      columnDefs: [{
          targets: pctargetsm,
          className: "dt-body-right",
          render: function(data,type,row){
              return parseFloat(data).toFixed(2);
          }
        }, {
          targets: [targetColm - 1],
          className: "dt-body-right",
          render: function(data,type,row){
              return parseFloat(data).toFixed(2) + '%';
          }
        }, {
          targets: [targetColm, targetColm+1, targetColm+2],
          className: "dt-body-center"
      }],
      buttons: [
          'copy', 'pdfHtml5','csvHtml5', 'colvis'
      ],
      colReorder: {order:colOrderm},
      select: true
  });

  $('#pcTableMFI').on('mouseover', 'tr', function() {
    var data = pcTablem.row(this).data();
    if (data != undefined) {
      var line = parseInt(data.idx);
      pcAppMFI.selectedLines = [line];
      updateParallelForegroundidx();
    }
  });
  $('#pcTableMFI').on('mouseleave', 'tr', function() {
    pcAppMFI.selectedLines = [];
    for (var i = 0, j = pcAppMFI.lines.length; i < j; i++){
      pcAppMFI.selectedLines.push(pcAppMFI.lines[i]);
    }
    updateParallelForegroundidx();
  });
};
/*
* Update Parallel Foreground
*/
var updateParallelForegroundidx = function() {
  pcAppMFI.foreground[0].map(function(d) {
    var ln = parseInt(d['__data__']['idx']);

    if ($.inArray(ln,pcAppMFI.selectedLines) < 0){
      d.style.display = "none";
    } else {
      d.style.display = null;
    }
  });
};
/*
 * Display The Main Plot
*/
var displayParallelPlotm = function() {
  var margin = {top: 30, right: 10, bottom: 10, left: 10},
      h = $("#chartDivPCm").height() * 0.6,
      w = $("#plotDivPCm").width(),
      y = {},
      dragging = {},
      width = w - margin.left - margin.right,
      height = h - margin.top - margin.bottom;

  $("#plotDivPCm").empty();
  $("#plotDivPCm").height(h);

  var svg = d3.select("#plotDivPCm").append("svg")
      .attr("width", w)
      .attr("height", h)
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
  var line = d3.svg.line();
  var axis = d3.svg.axis().orient("left").ticks(8);

  // Use this to scale line width to percentage population
  var pd = d3.extent(pcAppMFI.origData, function(p) { return +p['Percentage']; });
  var popScale = d3.scale.linear().range([1,5]).domain(pd);

  var dimensions = d3.keys(pcAppMFI.flowData[0]).filter(function(d) {
    return (y[d] = d3.scale.linear()
      .domain(d3.extent(pcAppMFI.flowData,function(p) { return +p[d]; }))
      .range([height, 0]));
  });
  x.domain(dimensions);

  function path(d) {
    return line(dimensions.map(function(p) { return [x(p), y[p](d[p])]; }));
  }
  function position(d) {
    var v = dragging[d];
    return v == null ? x(d) : v;
  }
  function transition(g) {
    return g.transition().duration(500);
  }

  function brush() {
    var actives = dimensions.filter(function(p) { return !y[p].brush.empty(); });
    var extents = actives.map(function(p) { return y[p].brush.extent(); });
    var indices = pcAppMFI.origData.filter(function(d) {
      var line = parseInt(d.idx)
      var tf = actives.every(function(p,i) {
        return extents[i][0] <= pcAppMFI.flowData[line][p] &&
                    pcAppMFI.flowData[line][p] <= extents[i][1];
      });
      if (tf) {
        return line.toString();
      }
    });

    pcAppMFI.selectedLines = indices.map(function(d) {
      return parseInt(d.idx);
    });
    pcAppMFI.lines = indices.map(function(d) {
      return parseInt(d.idx);
    });
    updateParallelForegroundidx();
    updatePopTableidx();
    updateSmpTableidx();
    displayTableGridm();
  }

  // Display paths in light gray color, to use as reference
  pcAppMFI.background = svg.append("g")
      .attr("class", "background")
    .selectAll("path")
      .data(pcAppMFI.flowData)
    .enter().append("path")
      .attr("d", path);

  // Add foreground lines for focus, color by population.
  pcAppMFI.foreground = svg.append("g")
      .attr("class", "foreground")
    .selectAll("path")
      .data(pcAppMFI.origData)
    .enter().append("path")
      .attr("d", path)
      .attr("stroke",function(d){
        var pop = pcAppMFI.populations.indexOf(d.Population) + 1;
        return color_palette[0][pop][0]; })
    //.attr("stroke-width", 1);
      // Use this if you want to scale the lines based on
      // population percentage
      .attr("stroke-width", function(d) {
        var pop = pcAppMFI.populations.indexOf(d.Population);
        var w = popScale(pcAppMFI.origData[pop]['Percentage']);
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
          pcAppMFI.background.attr("visibility", "hidden"); })
      .on("drag", function(d) {
          dragging[d] = Math.min(width, Math.max(0, d3.event.x));
          pcAppMFI.foreground.attr("d", path);
          dimensions.sort(function(a, b) { return position(a) - position(b); });
          x.domain(dimensions);
          g.attr("transform", function(d) { return "translate(" + position(d) + ")"; }); })
      .on("dragend", function(d) {
          delete dragging[d];
          transition(d3.select(this)).attr("transform", "translate(" + x(d) + ")");
          transition(pcAppMFI.foreground).attr("d", path);
          pcAppMFI.background
              .attr("d", path)
            .transition()
              .delay(500)
              .duration(0)
              .attr("visibility", null);
      }));

  // Add an axis and title.
  g.append("g")
      .attr("class", "axis")
      .each(function(d) { d3.select(this).call(axis.scale(y[d])); });
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
  $('#PCmline_opacity').on('change', (function() {
    var val = $(this).val();
    $('#plotDivPCm .foreground path').css('stroke-opacity', val.toString());
    $('#pcm_opacity').html((Math.round(val*10000)/100) + "%");
  }));
};

/*
 * Retrieve the data, then call display functions
*/
var displayParallelCoordinatesMFI = function() {
  var inputFile = "./csAllMFIs.tsv";
  d3.tsv(inputFile, function(error, data) {
    var allPops = 0,
        allSamples = 0;
    if (error) {
      alert("Problem Retrieving Data");
      return;
    }
    pcAppMFI.origData = $.extend(true,[],data);
    pcAppMFI.headers = Object.keys(pcAppMFI.origData[0]);
    pcAppMFI.headers.push("EditedPopName");
    pcAppMFI.origData.forEach(function(d,idx) {
      d.idx = idx;
      d.EditedPopName = d.Population;
      d.SmpName = d.SampleName;
      pcAppMFI.selectedLines.push(idx);
      pcAppMFI.lines.push(idx);
      if (!pcAppMFI.populations.includes(d.Population)){
        pcAppMFI.populations.push(d.Population);
      }
      if (!pcAppMFI.samples.includes(d.SmpName)){
        pcAppMFI.samples.push(d.SmpName);
      }
    });
    pcAppMFI.populations = pcAppMFI.populations.sort(function(a, b){return a-b});
    pcAppMFI.allLines = pcAppMFI.origData.length;

    allPops = pcAppMFI.populations.length;
    allSamples = pcAppMFI.samples.length;
    for (var i = 0; i < allPops; i++) {
      pcAppMFI.selectedPopulations.push(i);
    }
    for (var i = 0; i < allSamples; i++) {
      pcAppMFI.selectedSamples.push(i);
    }
    /*
     * For the plot use only the MFI information
     * for each populations. Store in flowData
    */
    pcAppMFI.flowData = $.extend(true,[],data);
    pcAppMFI.flowData.forEach(function(d) {
      delete d['Population'];
      delete d['SampleName'];
      delete d['Percentage'];
    });

    displayPopTablem();
    displaySmpTablem();
    displayTableGridm();
    displayParallelPlotm();

    $("#resetDisplayMFIpop").on("click",function() {
      var opcty = ".8";
      for (var i = 0; i < allPops; i++) {
        pcAppMFI.selectedPopulations.push(i);
      }
      for (var i = 0; i < allSamples; i++) {
        pcAppMFI.selectedSamples.push(i);
      }
      for (var i = 0; i < pcAppMFI.allLines; i++) {
        pcAppMFI.selectedLines.push(i);
        pcAppMFI.lines.push(i);
      }

      $("#popSelectAllPCm").prop('checked',true);
      $(".popSelectPCm").prop("checked",true);
      $("#smpSelectAllPCm").prop('checked',true);
      $(".smpSelectPCm").prop("checked",true);

      $('#plotDivPCm .foreground path').css('stroke-opacity', opcty);
      $('#pcm_opacity').html("80%");
      $('#PCmline_opacity').val(0.8);

      displayPopTablem();
      displaySmpTablem();
      displayTableGridm();
      displayParallelPlotm();
    });

    $(window).on('resize',function() {
      waitForFinalEvent(function() {
        displayAllm();
      }, 500, "resizePCm");
    });
  });
};
