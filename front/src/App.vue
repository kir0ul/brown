<template>
  <v-app>
    <v-tabs centered color="#3a1e1a" dark icons-and-text>
      <v-tabs-slider color="#c00404"></v-tabs-slider>

      <v-tab href="#tab-1">
        Articles
        <v-icon>description</v-icon>
      </v-tab>

      <v-tab href="#tab-2">
        Search
        <v-icon>search</v-icon>
      </v-tab>

      <v-tab href="#tab-3">
        Visualization
        <v-icon>insert_chart_outlined</v-icon>
      </v-tab>

      <v-tab-item v-for="i in 3" :value="'tab-' + i" :key="i">
        <v-content v-if="i === 1"> <DataTable /> </v-content>
      </v-tab-item>
    </v-tabs>
  </v-app>
</template>

<script>
import DataTable from "./components/DataTable";

export default {
  data() {
    return {};
  },
  components: {
    DataTable
  }
};
</script>
