// Copyright (c) 2016 Northrop Grumman.
// All rights reserved.

var processScatterDataMFI2D = function() {
  var xData = [],
      yData = [],
      popData = [],
      col1 = [],
      col2 = [],
      pop = scatterDataMFI['poplist'];

  col1 = scatterDataMFI['data'].map(function(value,index) {
    return value[scatterDataMFI['m1']];
  });
  col2 = scatterDataMFI['data'].map(function(value,index) {
    return value[scatterDataMFI['m2']];
  });
  for (var i = 0, j=col1.length; i < j; i++) {
    if (scatterDataMFI['selectedPopulations'].indexOf(pop[i]) >= 0) {
      xData.push(col1[i]);
      yData.push(col2[i]);
      popData.push(pop[i]);
    }
  }
  scatterDataMFI['popColors'] = popData.map(function(value,index) {
    return color_palette[0][value][0];
  });
  scatterDataMFI['xData'] = xData;
  scatterDataMFI['yData'] = yData;
  scatterDataMFI['popData'] = popData;
  return scatterDataMFI;
};

var processScatterData3DMFI = function() {
  var xData = [],
      yData = [],
      zData = [],
      col1 = [],
      col2 = [],
      col3 = [],
      pop = [],
      min = Number.MAX_VALUE,
      max = Number.MIN_VALUE,
      popData = [];

  min = d3.min(scatterData3DMFI['data'], function(array) {
    return d3.min(array);
  });
  max = d3.max(scatterData3DMFI['data'], function(array) {
    return d3.max(array);
  });
  scatterData3DMFI['min'] = 0;
  scatterData3DMFI['max'] = max;

  col1 = scatterData3DMFI['data'].map(function(value,index) {
    return value[scatterData3DMFI['m1']];
  });
  col2 = scatterData3DMFI['data'].map(function(value,index) {
    return value[scatterData3DMFI['m2']];
  });
  col3 = scatterData3DMFI['data'].map(function(value,index) {
    return value[scatterData3DMFI['m3']];
  });
  pop = scatterData3DMFI['poplist'];

  for (var i = 0, j = col1.length; i < j; i++) {
    if (scatterData3DMFI['selectedPopulations'].indexOf(pop[i]) >= 0) {
      xData.push(col1[i]);
      yData.push(col2[i]);
      zData.push(col3[i]);
      popData.push(pop[i]);
    }
  }

  scatterData3DMFI['popColors'] = popData.map(function(value,index) {
    return color_palette[0][value][0];
  });
  scatterData3DMFI['xData'] = xData;
  scatterData3DMFI['yData'] = yData;
  scatterData3DMFI['zData'] = zData;
  scatterData3DMFI['popData'] = popData;
  return scatterData3DMFI;
};
