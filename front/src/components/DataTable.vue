<template>
  <v-data-table
    :headers="headers"
    :items="itemData"
    hide-actions
    class="elevation-1"
  >
    <template slot="items" slot-scope="props">
      <td>{{ props.item.PMID }}</td>
      <td class="text-xs-left">{{ props.item.Title }}</td>
      <!-- <td class="text-xs-right">{{ props.item.Authors }}</td> -->
      <td class="text-xs-right">{{ props.item.PublicationYear }}</td>
    </template>
  </v-data-table>
</template>

<script>
import axios from "axios";
const url = "http://localhost:3000/articles";

export default {
  data() {
    return {
      headers: [
        {
          text: "PMID",
          align: "left",
          sortable: true,
          value: "PMID"
        },
        { text: "Title", align: "center", sortable: true, value: "Title" },
        // { text: "Authors", align: "center", sortable: true, value: "Authors" },
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
  mounted() {
    axios({
      method: "GET",
      url: url
    }).then(
      response => {
        this.itemData = response.data;
      },
      error => {
        console.error(error);
      }
    );
  }
};
</script>
