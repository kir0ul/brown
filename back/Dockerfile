FROM node:carbon-slim

WORKDIR /usr/src/app

COPY package*.json ./

RUN npm install

COPY . .

EXPOSE 3000

ENV NODE_ENV docker
CMD ["npm", "run", "start"]
