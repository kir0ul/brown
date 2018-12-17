<template>
  <v-container fluid>
    <v-layout row wrap align-center justify-center>
      <v-flex xs6 text-xs-center>
        <v-btn block v-on:click="createHorizontalBarChart">Draw</v-btn>
      </v-flex>
      <v-flex xs12> <div class="js-bar-container"></div> </v-flex>
    </v-layout>
  </v-container>
</template>

<script>
// import LineChart from "britecharts/dist/umd/line.min";
import BarChart from "britecharts/dist/umd/bar.min";
const d3Selection = require("d3-selection");

export default {
  data() {
    return {
      data: [
        { name: "Shiny", id: 1, quantity: 86, percentage: 5 },
        { name: "Blazing", id: 2, quantity: 300, percentage: 18 },
        { name: "Dazzling", id: 3, quantity: 276, percentage: 16 },
        { name: "Radiant", id: 4, quantity: 195, percentage: 11 },
        { name: "Sparkling", id: 5, quantity: 36, percentage: 2 },
        { name: "Other", id: 0, quantity: 814, percentage: 48 }
      ]
    };
  },
  methods: {
    createHorizontalBarChart: function() {
      let barChart = new BarChart(),
        margin = {
          left: 120,
          right: 20,
          top: 20,
          bottom: 30
        },
        barContainer = d3Selection.select(".js-bar-container"),
        containerWidth = barContainer.node()
          ? barContainer.node().getBoundingClientRect().width
          : false;
      barChart
        .margin(margin)
        .width(containerWidth)
        .valueLabel("percentage")
        .height(300);
      // .colorSchema(britecharts.colors.colorSchemas.britechartsColorSchema)

      barContainer.datum(this.data).call(barChart);
    }
  }
};
</script>

<link
  type="text/css"
  rel="stylesheet"
  href="node_modules/britecharts/dist/css/britecharts.min.css"
/>
