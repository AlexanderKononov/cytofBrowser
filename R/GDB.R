
#### Create api for neo4j session
get_neo_api <- function(user = "neo4j", password = "password"){
  neo_api <- neo4j_api$new(url = "http://localhost:7474", user = user, password = password)
}
