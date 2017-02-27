// Copyright (c) 2016 Northrop Grumman.
// All rights reserved.

var displayMarkerTable = function(plotconfig){
  var nbm = plotconfig.mrkrNames.length + 1;
  $(plotconfig.mtable).empty();
  plotconfig.allMarkers.map(function(v) {
    $(plotconfig.mtable).append('<tr>'
        + '<td><span style="background-color:rgba(0,0,0,' + (v + 1 )/ nbm + ')'
        + '">&nbsp;&nbsp;&nbsp;</span></td>'
        + '<td title="' + plotconfig.mrkrNames[v] + '">'
        + plotconfig.mrkrNames[v] + '</td>'
        + '<td align="center"><input type="checkbox" checked class='
        + plotconfig.mrkrSelect + ' value=' + v + '/></td></tr>');
  });
  if (nbm > 5) {
    $(plotconfig.mrkrSelectAll).prop('checked', false);
    $(plotconfig.mrkrSelectAll).prop('disabled', true);
    $('#markerWarning').show();
    $(plotconfig.mrkrSelectj).each(function() {
      var selectedMrkr = parseInt(this.value);
      if (selectedMrkr > 4){
        this.checked = false;
        this.disabled = true;
      } else {
        this.checked = true;
      }
    });
  }

  $(plotconfig.mrkrSelectAll).click(function() {
    var checkAll = $(plotconfig.mrkrSelectAll).prop('checked');
    if (checkAll) {
      $(plotconfig.mrkrSelectj).prop("checked", true);
    } else {
      $(plotconfig.mrkrSelectj).prop("checked", false);
    }
    updateCSplots(plotconfig);
  });

  $(plotconfig.mrkrSelectj).click(function() {
    if (nbm < 6){
      if ($(plotconfig.mrkrSelectj).length == $(plotconfig.mrkrSelectCheck).length) {
        $(plotconfig.mrkrSelectAll).prop("checked",true);
      } else {
        $(plotconfig.mrkrSelectAll).prop("checked",false);
      }
    } else {
      var nbSelected = 0;
      $(plotconfig.mrkrSelectj).each(function() {
        if (this.checked) {nbSelected++}
      });
      if (nbSelected < 5) {
        $(plotconfig.mrkrSelectj).prop('disabled', false);
      } else {
        $(plotconfig.mrkrSelectj).each(function() {
          if (!this.checked) {
            this.disabled = true;
          }
        });
      }
    }
    updateCSplots(plotconfig);
  });
};

var updateBoxplot = function(plotconfig){
  var margin = {top: 30, right: 10, bottom: 50, left: 60},
      h = 0,
      w = 0,
      width = 0,
      height = 0,
      labels = false, // show the text labels beside individual boxplots?
      mfi_option = false,
      min = Infinity,
      max = -Infinity,
      checkLabels = $(plotconfig.displayvalues).prop("checked"),
      checkMFI = $(plotconfig.displayMFI).prop("checked"),
      dataToPlot = [],
      tmp = [],
      nbm = plotconfig.mrkrNames.length + 1,
      maxRange = 0,
      minRange = 0,
      domainx = [],
      domainx1 = [];

  $(plotconfig.plotdivj).empty();
  h = $(window).height() - 200;
  $(plotconfig.plotdivj).height(h);
  w = $(plotconfig.plotdivj).width();
  width = w - margin.left - margin.right;
  height = h - margin.top - margin.bottom;

  var svg = d3.select(plotconfig.plotdivj).append("svg")
      .attr("width", w)
      .attr("height", h)
      .attr("class", "box")
    .append("g")
      .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

  if (checkLabels) {
      labels = true;
  };
  if (checkMFI) {
      mfi_option = true;
  };
  /* Get the data in proper shape to feed to the boxplot function
  want [Object, Object, ..., Object] where Object:
    'population': pop
    'Data' : [Object, Object, ..., Object] where Object:
        'Marker' : marker name
        'outliers' : outliers
  */

  if (plotconfig.view == 'p') {
    plotconfig.selectedPopulations.forEach(function(p) {
      tmpPlot = [];
      plotconfig.selectedMarkers.forEach(function(m) {
        var markernm = plotconfig.mrkrNames[m],
            qtmp = [
                +plotconfig.csdata.q1[markernm][p],
                +plotconfig.csdata.q2[markernm][p],
                +plotconfig.csdata.q3[markernm][p]
            ],
            wtmp = [
                +plotconfig.csdata.lower[markernm][p],
                +plotconfig.csdata.upper[markernm][p]
            ];
        tmp = [];
        // Get min and max while we're here
        plotconfig.csdata.outliers[markernm][p].forEach(function(outl) {
          tmp.push(+outl);
          if (+outl > max) {max = +outl};
          if (+outl < min) {min = +outl};
        });
        if (+plotconfig.csdata.upper[markernm][p] > max) {
          max = +plotconfig.csdata.upper[markernm][p];
        };
        if (+plotconfig.csdata.lower[markernm][p] < min) {
          min = +plotconfig.csdata.lower[markernm][p];
        };
        tmpPlot.push({
            marker: markernm,
            outliers: tmp,
            quartiles: qtmp,
            whiskers: wtmp,
            config: [m,p, nbm],
            mfi: +plotconfig.csdata.mfi[markernm][p]
        });
      });
      dataToPlot.push({population:p, popdata: tmpPlot});
    });
  } else {
    plotconfig.selectedMarkers.forEach(function(m) {
      var markernm = plotconfig.mrkrNames[m];
      tmpPlot = [];
      plotconfig.selectedPopulations.forEach(function(p) {
        var qtmp = [
                +plotconfig.csdata.q1[markernm][p],
                +plotconfig.csdata.q2[markernm][p],
                +plotconfig.csdata.q3[markernm][p]
            ],
           wtmp = [
                +plotconfig.csdata.lower[markernm][p],
                +plotconfig.csdata.upper[markernm][p]
            ];
        tmp = [];
        // Get min and max while we're here
        plotconfig.csdata.outliers[markernm][p].forEach(function(outl) {
          tmp.push(+outl);
          if (+outl > max) {max = +outl};
          if (+outl < min) {min = +outl};
        });
        if (+plotconfig.csdata.upper[markernm][p] > max) {
          max = +plotconfig.csdata.upper[markernm][p];
        };
        if (+plotconfig.csdata.lower[markernm][p] < min) {
          min = +plotconfig.csdata.lower[markernm][p];
        };
        tmpPlot.push({
            population:p,
            outliers: tmp,
            quartiles: qtmp,
            whiskers: wtmp,
            config: [m,p, nbm],
            mfi: +plotconfig.csdata.mfi[markernm][p]
        });
      });
      dataToPlot.push({marker: markernm, popdata: tmpPlot});
    });
  };
  maxRange = max + 30;
  minRange = min - 30;

  if (plotconfig.view == 'p') {
    domainx = plotconfig.selectedPopulations;
    domainx1 = plotconfig.selectedMarkers.map(function(d){
      return plotconfig.mrkrNames[d];
    });
  } else {
    domainx1 = plotconfig.selectedPopulations;
    domainx = plotconfig.selectedMarkers.map(function(d){
      return plotconfig.mrkrNames[d];
    });
  }
  // axes
  var xScale = d3.scale.ordinal()
      .domain(domainx)
      .rangeRoundBands([0 , width], 0.2, 0.02);

  var x1Scale = d3.scale.ordinal()
      .domain(domainx1)
      .rangeRoundBands([0, xScale.rangeBand()], 0.1);

  var xAxis = d3.svg.axis()
      .scale(xScale)
      .orient("bottom");

  // the y-axis
  var yScale = d3.scale.linear()
      .domain([minRange, maxRange])
      .range([height + margin.top, 0 + margin.top]);

  var yAxis = d3.svg.axis()
      .scale(yScale)
      .orient("left")
      .tickFormat(d3.format("d"));

  svg.append("g")
      .attr("class", "x axisbp")
      .attr("transform", "translate(0," + (height + margin.top) + ")")
      .call(xAxis);

  svg.append("g")
      .attr("class", "y axisbp")
      .call(yAxis)
    .append("text")
      .attr("class", "ylabel")
      .attr("transform", "rotate(-90)")
      .attr("y", 0 - margin.left)
      .attr("x", 0 - (height / 2))
      .attr("dy", "1em")
      .style("text-anchor", "middle")
      .text("MFI values");

  var boxplot = d3.box()
      .width(x1Scale.rangeBand())
      .height(height + margin.top)
      .domain([minRange, maxRange])
      .showLabels(labels)
      .showMFI(mfi_option);

  if (plotconfig.view == 'p'){
    var group = svg.selectAll(".groups")
        .data(dataToPlot)
      .enter().append("g")
        .attr("class", "group")
        .attr("transform", function(d) {
            return "translate(" + xScale(d.population) + ",0)";
        });

    group.selectAll(".box")
        .data(function(d) { return d.popdata; })
      .enter().append("g")
        .attr("transform", function(d) {return "translate(" + x1Scale(d.marker) + ",0)"; })
        .call(boxplot);
  } else {
    var group = svg.selectAll(".groups")
        .data(dataToPlot)
      .enter().append("g")
        .attr("class", "group")
        .attr("transform", function(d) {
            return "translate(" + xScale(d.marker) + ",0)";
        });

      group.selectAll(".box")
          .data(function(d) { return d.popdata; })
        .enter().append("g")
          .attr("transform", function(d) { return "translate(" + x1Scale(d.population) + ",0)"; })
          .call(boxplot);
  }
};

(function() {
  // Inspired by http://informationandvisualization.de/blog/box-plot
  // Modified to fit our data structure.
  d3.box = function() {
    var width = 1,
        height = 1,
        duration = 0,
        domain = null,
        value = Number,
        showLabels = true, // whether or not to show text labels
        numBars = 4,
        curBar = 1,
        showMFI = true, // display MFI ?
        tickFormat = null,
        margin = {top: 30, right: 10, bottom: 50, left: 60};

    // For each small multiple…
    function box(g) {
      g.each(function(data, i) {
        var d = data.outliers.sort(d3.ascending),
            g = d3.select(this),
            n = d.length,
            min = Infinity,
            max = -Infinity;
        if (n > 0){
            min = d[0],
            max = d[n - 1];
        }
        // Readjust min and max with upper and lower values
        if (data.whiskers[0] < min) {min = data.whiskers[0]}
        if (data.whiskers[1] > max) {max = data.whiskers[1]}
        // Compute quartiles. Must return exactly 3 elements.
        var quartileData = data.quartiles;
        // Compute whiskers. Must return exactly 2 elements, or null.
        var whiskerData = data.whiskers;
        // Compute outliers. here all data in d is an outlier.
        // We compute the outliers as indices, so that we can join across transitions!
        var outlierIndices = d3.range(n);
        var mfiData = data.mfi;
        // this is the scale for ONE SET of values
        // Compute the new x-scale.
        var x1 = d3.scale.linear()
            .domain(domain && domain.call(this, d, i) || [min, max])
            .range([height , 0 + margin.top ]);
        // Retrieve the old x-scale, if this is an update.
        var x0 = this.__chart__ || d3.scale.linear()
            .domain([0, Infinity])
            .range(x1.range());

        // Stash the new scale.
        this.__chart__ = x1;
// Note: the box, median, and box tick elements are fixed in number,
// so we only have to handle enter and update. In contrast, the outliers
// and other elements are variable, so we need to exit them! Variable
// elements also fade in and out.
        // Update center line: the vertical line spanning the whiskers.
        var center = g.selectAll("line.center")
            .data(whiskerData ? [whiskerData] : []);

        //vertical line
        center.enter().insert("line", "rect")
            .attr("class", "center")
            .attr("x1", width / 2)
            .attr("y1", function(d) { return x0(d[0]); })
            .attr("x2", width / 2)
            .attr("y2", function(d) { return x0(d[1]); })
            .style("opacity", 1e-6)
            .style("stroke", function(d) { return color_palette[0][data.config[1]][3]; })
          .transition()
            .duration(duration)
            .style("opacity", 1)
            .attr("y1", function(d) { return x1(d[0]); })
            .attr("y2", function(d) { return x1(d[1]); });

        center.transition()
            .duration(duration)
            .style("opacity", 1)
            .attr("y1", function(d) { return x1(d[0]); })
            .attr("y2", function(d) { return x1(d[1]); });

        center.exit().transition()
            .duration(duration)
            .style("opacity", 1e-6)
            .attr("y1", function(d) { return x1(d[0]); })
            .attr("y2", function(d) { return x1(d[1]); })
            .remove();

        // Update innerquartile box.
        var box = g.selectAll("rect.box")
            .data([quartileData]);

        box.enter().append("rect")
            .attr("class", "box")
            .style("fill", function(d) {
                var nbm = data.config[2],
                    pop = data.config[1],
                    mrkr = data.config[0];
                var color = color_palette[0][pop][1] + (mrkr + 1 )/ nbm + ')';
                return color; })
            .style("stroke", function(d) { return color_palette[0][data.config[1]][3]; })
            .attr("x", 0)
            .attr("y", function(d) { return x0(d[2]); })
            .attr("width", width)
            .attr("height", function(d) { return x0(d[0]) - x0(d[2]); })
          .transition()
            .duration(duration)
            .attr("y", function(d) { return x1(d[2]); })
            .attr("height", function(d) { return x1(d[0]) - x1(d[2]); });

        box.transition()
            .duration(duration)
            .attr("y", function(d) { return x1(d[2]); })
            .attr("height", function(d) { return x1(d[0]) - x1(d[2]); });

        // Update median line.
        var medianLine = g.selectAll("line.median")
            .data([quartileData[1]]);

        medianLine.enter().append("line")
            .attr("class", "median")
            .attr("x1", 0)
            .attr("y1", x0)
            .attr("x2", width)
            .attr("y2", x0)
            .style("stroke", function(d) { return color_palette[0][data.config[1]][3]; })
          .transition()
            .duration(duration)
            .attr("y1", x1)
            .attr("y2", x1);

        medianLine.transition()
            .duration(duration)
            .attr("y1", x1)
            .attr("y2", x1);

        // Update MFI line.
        var MFILine = g.selectAll("line.mfi")
            .data([mfiData]);
        if (showMFI == true) {
            MFILine.enter().append("line")
                .attr("class", "mfi")
                .style("stroke", function(d){ return color_palette[0][data.config[1]][2]; })
                .attr("x1", 0)
                .attr("y1", x0)
                .attr("x2", width)
                .attr("y2", x0)
              .transition()
                .duration(duration)
                .attr("y1", x1)
                .attr("y2", x1);

            MFILine.transition()
                .duration(duration)
                .attr("y1", x1)
                .attr("y2", x1);
        }

        // Update whiskers.
        var whisker = g.selectAll("line.whisker")
            .data(whiskerData || []);

        whisker.enter().insert("line", "circle, text")
            .attr("class", "whisker")
            .attr("x1", 0)
            .attr("y1", x0)
            .attr("x2", 0 + width)
            .attr("y2", x0)
            .style("opacity", 1e-6)
            .style("stroke", function(d) { return color_palette[0][data.config[1]][3]; })
          .transition()
            .duration(duration)
            .attr("y1", x1)
            .attr("y2", x1)
            .style("opacity", 1);

        whisker.transition()
            .duration(duration)
            .attr("y1", x1)
            .attr("y2", x1)
            .style("opacity", 1);

        whisker.exit().transition()
            .duration(duration)
            .attr("y1", x1)
            .attr("y2", x1)
            .style("opacity", 1e-6)
            .remove();

        // Update outliers.
        var outlier = g.selectAll("circle.outlier")
            .data(outlierIndices, Number);

        outlier.enter().insert("circle", "text")
            .attr("class", "outlier")
            .attr("r", 3)
            .attr("cx", function(d){
                  return Math.floor(Math.random() * width);
                })
            .attr("cy", function(i) { return x0(d[i]); })
            .style("opacity", 1e-6)
            .style("fill", function(d) {
                    var nbm = data.config[2],
                        pop = data.config[1],
                        mrkr = data.config[0];
                    var color = color_palette[0][pop][1] + (mrkr + 1 )/ nbm + ')';
                    return color; })
            .style("stroke", function(d) { return color_palette[0][data.config[1]][3]; })
          .transition()
            .duration(duration)
            .attr("cy", function(i) { return x1(d[i]); })
            .style("opacity", 1);

        outlier.transition()
            .duration(duration)
            .attr("cy", function(i) { return x1(d[i]); })
            .style("opacity", 1);

        outlier.exit().transition()
            .duration(duration)
            .attr("cy", function(i) { return x1(d[i]); })
            .style("opacity", 1e-6)
            .remove();

        // Compute the tick format.
        var format = tickFormat || x1.tickFormat(8);

        // Update box ticks.
        var boxTick = g.selectAll("text.box")
            .data(quartileData);
        if(showLabels == true) {
          boxTick.enter().append("text")
              .attr("class", "box")
              .attr("dy", ".3em")
              .attr("dx", function(d, i) { return i & 1 ? 6 : -6 })
              .attr("x", function(d, i) { return i & 1 ?  + width : 0 })
              .attr("y", x0)
              .attr("text-anchor", function(d, i) { return i & 1 ? "start" : "end"; })
              .text(format)
            .transition()
              .duration(duration)
              .attr("y", x1);
        }

        boxTick.transition()
            .duration(duration)
            .text(format)
            .attr("y", x1);

        // Update whisker ticks. These are handled separately from the box
        // ticks because they may or may not exist, and we want don't want
        // to join box ticks pre-transition with whisker ticks post-.
        var whiskerTick = g.selectAll("text.whisker")
            .data(whiskerData || []);
        if(showLabels == true) {
          whiskerTick.enter().append("text")
              .attr("class", "whisker")
              .attr("dy", ".3em")
              .attr("dx", 6)
              .attr("x", width)
              .attr("y", x0)
              .text(format)
              .style("opacity", 1e-6)
            .transition()
              .duration(duration)
              .attr("y", x1)
              .style("opacity", 1);
        }
        whiskerTick.transition()
            .duration(duration)
            .text(format)
            .attr("y", x1)
            .style("opacity", 1);

        whiskerTick.exit().transition()
            .duration(duration)
            .attr("y", x1)
            .style("opacity", 1e-6)
            .remove();
      });
      d3.timer.flush();
    }

    box.width = function(x) {
      if (!arguments.length) return width;
      width = x;
      return box;
    };
    box.height = function(x) {
      if (!arguments.length) return height;
      height = x;
      return box;
    };
    box.tickFormat = function(x) {
      if (!arguments.length) return tickFormat;
      tickFormat = x;
      return box;
    };
    box.duration = function(x) {
      if (!arguments.length) return duration;
      duration = x;
      return box;
    };
    box.domain = function(x) {
      if (!arguments.length) return domain;
      domain = x == null ? x : d3.functor(x);
      return box;
    };
    box.value = function(x) {
      if (!arguments.length) return value;
      value = x;
      return box;
    };
    box.showLabels = function(x) {
      if (!arguments.length) return showLabels;
      showLabels = x;
      return box;
    };
    box.showMFI = function(x) {
      if (!arguments.length) return showMFI;
      showMFI = x;
      return box;
    };
    return box;
  };
})();
