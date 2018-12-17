<template>
  <v-container fluid>
    <v-layout row wrap align-center justify-center>
      <v-flex xs4>
        <v-text-field
          ref="SearchContent"
          label="Search author"
          v-model="DefaultText"
          clearable
          v-on:keyup.enter="searchAuthor"
        ></v-text-field>
      </v-flex>
      <v-flex xs1>
        <v-btn flat icon v-on:click="searchAuthor">
          <v-icon>search</v-icon>
        </v-btn>
      </v-flex>
    </v-layout>

    <v-layout row wrap align-center>
      <v-flex xs12>
        <v-data-table
          :headers="headers"
          :items="itemData"
          hide-actions
          class="elevation-1"
        >
          <template slot="items" slot-scope="props">
            <td>{{ props.item.PMID }}</td>
            <td class="text-xs-left">{{ props.item.Title }}</td>
            <td class="text-xs-right">{{ props.item.LastName }}</td>
            <td class="text-xs-right">{{ props.item.PublicationYear }}</td>
          </template>
        </v-data-table>
      </v-flex>
    </v-layout>
  </v-container>
</template>

<script>
import axios from "axios";
const url = "http://localhost:3000/author/";

export default {
  data() {
    return {
      DefaultText: "",
      headers: [
        {
          text: "PMID",
          align: "left",
          sortable: true,
          value: "PMID"
        },
        { text: "Title", align: "center", sortable: true, value: "Title" },
        { text: "Author", align: "center", sortable: true, value: "LastName" },
        {
          text: "Publication year",
          align: "center",
          sortable: true,
          value: "PublicationYear"
        }
      ],
      itemData: []
    };
  },
  methods: {
    searchAuthor: function() {
      let SearchedAuthor = this.$refs.SearchContent.value;
      if (SearchedAuthor) {
        axios({
          method: "GET",
          url: url + SearchedAuthor
        }).then(
          response => {
            this.itemData = response.data;
          },
          error => {
            console.error(error);
          }
        );
      }
    }
  }
};
</script>
