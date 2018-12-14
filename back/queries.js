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
  "article",
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
const Authors = sequelize.define(
  "author",
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

////////// Queries //////////
const getArticles = (request, response) =>
  Article.findAll({
    order: [
      sequelize.literal('"PublicationYear" DESC NULLS LAST'),
      sequelize.literal('"PublicationMonth" DESC NULLS LAST'),
      sequelize.literal('"PublicationDay" DESC NULLS LAST')
    ],
    limit: 500
  }).then(results => {
    response.status(200).json(results);
  });

const getAuthors = (request, response) =>
  Authors.findAll({ limit: 20 }).then(results => {
    response.status(200).json(results);
  });

module.exports = {
  getArticles,
  getAuthors
};
