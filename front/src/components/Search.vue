<template>
  <v-form>
    <v-container fluid>
      <v-layout row wrap align-center>
        <v-flex xs4 offset-xs4>
          <v-text-field
            v-on:keyup.enter="searchAuthor()"
            label="Search author"
            append-icon="search"
            v-model="DefaultText"
            clearable
          ></v-text-field>
        </v-flex>

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
              <td class="text-xs-right">{{ props.item.Authors }}</td>
              <td class="text-xs-right">{{ props.item.PublicationYear }}</td>
            </template>
          </v-data-table>
        </v-flex>
      </v-layout>
    </v-container>
  </v-form>
</template>

<script>
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
        { text: "Authors", align: "center", sortable: true, value: "Authors" },
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
      console.log(this.value);
    }
  }
};
</script>
