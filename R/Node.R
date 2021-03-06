#' Obtain a node in mother class object
#'
#' @description Use a path (numeric vector) to obtain a subgroup of a structure (mother class object).
#'
#' @param path the path of the node (numeric vector).
#' @param structure a mother class object (S4).
#'
#' @details Every node of a mother object (structure) can be identified with a numeric vector that indicates
#' the path used from the root to the node. The vector is the 'path' argument and is used to find specific
#' nodes of a given structure. For a complete explanation, we refer to Cossette et al. (2017).
#'
#' @examples
#' # We directly give the path of the desired node.
#' Node(c(0,2,2), LOG(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                               LOG(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                               GAMMA(1/30, c(3,4), NULL))))))
#'
#' # Here we provide the path with the GeneticCodes function of this package.
#' structure <- LOG(0.5, NULL, list(GAMMA(1/30, c(5,6), NULL),
#'                                      LOG(0.1, NULL, list(GAMMA(1/30, c(1,2), NULL),
#'                                      GAMMA(1/30, c(3,4), NULL)))))
#' Node(GeneticCodes(structure)[[3]], structure)
#'
#' @return Either a child or mother class object.
#'
#' @author Simon-Pierre Gadoury
#'
#' @export

Node <- function(path, structure)
{
  if (length(path) == 1)
    structure
  else
  {
    if (path[2] == 0)
      structure
    else
    {
      struc <- structure@structure[[path[2]]]
      Node(path[-2], struc)
    }
  }
}
