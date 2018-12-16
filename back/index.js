const express = require("express");
const bodyParser = require("body-parser");
const db = require("./queries");

// Run server
const app = express();
const port = 3000;
app.use(bodyParser.json());
app.use(
  bodyParser.urlencoded({
    extended: true
  })
);

// Enable CORS
app.use(function(req, res, next) {
  res.header("Access-Control-Allow-Origin", "*");
  res.header(
    "Access-Control-Allow-Headers",
    "Origin, X-Requested-With, Content-Type, Accept"
  );
  next();
});

// API root
app.get("/", (request, response) => {
  response.json({ info: "API root" });
});

// API routes
app.get("/articles", db.getArticles);
app.get("/author/:lastname", db.SearchAuthor);

// No route found
app.use(function(req, res) {
  res.status(404).send({ url: req.originalUrl + " not found" });
});

app.listen(port, () => {
  console.log(`App running on port ${port}.`);
});
