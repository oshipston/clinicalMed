var http = require('http');

var server = http.createServer(function(req, res) {
res.writeHead(200);
res.end('Hi everybody!');
});
server.listen(8080);

// nohup node http_server.js > output.log &  // starts server as background process

// Should consider use of 'npm install forever -g'
// capable of restarting closed instances of the server
// forever start http_server.js
// forever list
// forever stop http_server.js
