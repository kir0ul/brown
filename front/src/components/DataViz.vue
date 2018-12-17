<template>
  <v-container fluid>
    <v-layout row wrap align-center justify-center>
      <v-flex xs6 text-xs-center>
        <v-btn block v-on:click="createHorizontalBarChart">Draw</v-btn>
      </v-flex>
      <v-flex xs12> <div class="js-bar-container"></div> </v-flex>
      <v-flex xs6 text-xs-center>
        <h3 v-show="show">Number of publications by year</h3>
      </v-flex>
    </v-layout>
  </v-container>
</template>

<script>
// import LineChart from "britecharts/dist/umd/line.min";
import BarChart from "britecharts/dist/umd/bar.min";
import colors from "britecharts/dist/umd/colors.min";
const d3Selection = require("d3-selection");
import axios from "axios";
const url = "http://localhost:3000/dataviz";

export default {
  data() {
    return {
      show: false
    };
  },
  methods: {
    createHorizontalBarChart: function() {
      axios({
        method: "GET",
        url: url
      }).then(
        response => {
          var data = response.data;
          var dataset = [];
          for (let i = 0; i < data.length; i++) {
            dataset.push({
              val: Number(data[i].value),
              name: String(data[i].name)
            });
          }

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
            .colorSchema(colors.colorSchemas.britecharts)
            .valueLabel("val")
            .height(300)
            .yAxisLabel("Publications number")
            .xAxisLabel("Year");

          barContainer.datum(dataset).call(barChart);

          this.show = true;
        },
        error => {
          console.error(error);
        }
      );
    }
  }
};
</script>

<link
  type="text/css"
  rel="stylesheet"
  href="node_modules/britecharts/dist/css/britecharts.min.css"
/>
