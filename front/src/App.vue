<template>
  <v-app>
    <v-tabs centered color="#3a1e1a" dark icons-and-text>
      <v-tabs-slider color="#c00404"></v-tabs-slider>

      <v-tab href="#tab-1">
        Home
        <v-icon>home</v-icon>
      </v-tab>

      <v-tab href="#tab-2">
        Articles
        <v-icon>description</v-icon>
      </v-tab>

      <v-tab href="#tab-3">
        Search authors
        <v-icon>search</v-icon>
      </v-tab>

      <v-tab href="#tab-4">
        Visualization
        <v-icon>insert_chart_outlined</v-icon>
      </v-tab>

      <v-tab-item v-for="i in 4" :value="'tab-' + i" :key="i">
        <v-content v-if="i === 2"> <DataTable /> </v-content>
        <v-content v-if="i === 3"> <Search /> </v-content>
        <v-content v-if="i === 4"> <DataViz /> </v-content>
      </v-tab-item>
    </v-tabs>
  </v-app>
</template>

<script>
import DataTable from "./components/DataTable";
import Search from "./components/Search";
import DataViz from "./components/DataViz";

export default {
  data() {
    return {};
  },
  components: {
    DataTable,
    Search,
    DataViz
  }
};
</script>
