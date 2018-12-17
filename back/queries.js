const Sequelize = require("sequelize");
var config = require("config");

////////// Database config //////////
var dbName = config.get("database.name");
var dbUsername = config.get("database.username");
var dbPassword = config.get("database.password");
var dbHost = config.get("database.host");
var dbDialect = config.get("database.dialect");

const sequelize = new Sequelize(dbName, dbUsername, dbPassword, {
  host: dbHost,
  dialect: dbDialect,
  operatorsAliases: false,

  pool: {
    max: 5,
    min: 0,
    acquire: 30000,
    idle: 10000
  }
});

sequelize
  .authenticate()
  .then(() => {
    console.log("Connection has been established successfully.");
  })
  .catch(err => {
    console.error("Unable to connect to the database:", err);
  });

////////// Database model //////////
const Article = sequelize.define(
  "articles",
  {
    id: { type: Sequelize.INTEGER, primaryKey: true },
    PMID: Sequelize.STRING,
    Title: Sequelize.STRING,
    PublicationYear: Sequelize.INTEGER,
    PublicationMonth: Sequelize.INTEGER,
    PublicationDay: Sequelize.INTEGER
  },
  {
    // don't add the timestamp attributes (updatedAt, createdAt)
    timestamps: false
  }
);
const Author = sequelize.define(
  "authors",
  {
    id: { type: Sequelize.INTEGER, primaryKey: true },
    LastName: Sequelize.STRING,
    ForeName: Sequelize.STRING,
    Initials: Sequelize.INTEGER,
    ArticleId: Sequelize.INTEGER
  },
  {
    // don't add the timestamp attributes (updatedAt, createdAt)
    timestamps: false
  }
);
Article.hasMany(Author, { foreignKey: "ArticleId" });
Author.belongsTo(Article, { foreignKey: "id" });

/* const getArticles = (request, response) =>
 *   Article.findAll({
 *     order: [
 *       sequelize.literal('"PublicationYear" DESC NULLS LAST'),
 *       sequelize.literal('"PublicationMonth" DESC NULLS LAST'),
 *       sequelize.literal('"PublicationDay" DESC NULLS LAST')
 *     ],
 *     limit: 20
 *   })
 *     .then(results => {
 *       response.status(200).json(results);
 *     })
 *     .catch(err => {
 *       console.error(err);
 *       response.status(500).send("Error...");
 *     }); */
const getArticles = (request, response) =>
  sequelize
    .query(
      [
        "SELECT *",
        "FROM articles",
        "JOIN authors",
        'ON articles.id = authors."ArticleId"',
        "GROUP BY articles.id, authors.id",
        'ORDER BY "PublicationYear" DESC NULLS LAST,',
        '"PublicationMonth" DESC NULLS LAST,',
        '"PublicationDay" DESC NULLS LAST',
        "LIMIT 20"
      ].join("\n")
    )
    .then(results => {
      response.status(200).json(results);
    })
    .catch(err => {
      console.error(err);
      response.status(500).send("Error...");
    });

// Must do the query by hand because this dumb ORM is unable to do a JOIN with the IDs you want...
const SearchAuthor = (request, response) =>
  sequelize
    .query(
      [
        "SELECT *",
        "FROM authors",
        "JOIN articles",
        'ON articles.id = authors."ArticleId"',
        'WHERE LOWER("LastName") LIKE :queriedLastName',
        'ORDER BY "PublicationYear" DESC NULLS LAST,',
        '"PublicationMonth" DESC NULLS LAST,',
        '"PublicationDay" DESC NULLS LAST',
        "LIMIT 20"
      ].join("\n"),
      {
        replacements: {
          queriedLastName: "%" + request.params.lastname.toLowerCase() + "%"
        },
        type: sequelize.QueryTypes.SELECT
      }
    )
    .then(results => {
      response.status(200).json(results);
    })
    .catch(err => {
      console.error(err);
      response.status(500).send("Error...");
    });

const dataviz = (request, response) =>
  Article.findAll({
    attributes: [
      [sequelize.fn("COUNT", sequelize.col("id")), "articles_nb"],
      "PublicationYear",
      "PublicationMonth"
    ],
    group: ["PublicationYear", "PublicationMonth"],
    order: [
      sequelize.literal('"PublicationYear" ASC NULLS FIRST'),
      sequelize.literal('"PublicationMonth" ASC NULLS FIRST')
    ]
  })
    .then(results => {
      response.status(200).json(results);
    })
    .catch(err => {
      console.error(err);
      response.status(500).send("Error...");
    });

// Export queries
module.exports = {
  getArticles,
  SearchAuthor,
  dataviz
};
