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

      <v-tab-item v-for="i in 3" :id="'tab-' + i" :key="i">
        <v-card flat>
          <v-card-text>{{ text }}</v-card-text>
        </v-card>
      </v-tab-item>
    </v-tabs>
  </v-app>
</template>

<script>
export default {
  data() {
    return {
      text:
        "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat."
    };
  }
};
</script>
