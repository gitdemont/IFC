################################################################################
# This file is released under the GNU General Public License, Version 3, GPL-3 #
# Copyright (C) 2025 Yohann Demont                                             #
#                                                                              #
# It is part of IFC package, please cite:                                      #
# -IFC: An R Package for Imaging Flow Cytometry                                #
# -YEAR: 2020                                                                  #
# -COPYRIGHT HOLDERS: Yohann Demont, Gautier Stoll, Guido Kroemer,             #
#                     Jean-Pierre Marolleau, Loïc Garçon,                      #
#                     INSERM, UPD, CHU Amiens                                  #
#                                                                              #
# DISCLAIMER:                                                                  #
# -You are using this package on your own risk!                                #
# -We do not guarantee privacy nor confidentiality.                            #
# -This program is distributed in the hope that it will be useful, but WITHOUT #
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or        #
# FITNESS FOR A PARTICULAR PURPOSE. In no event shall the copyright holders or #
# contributors be liable for any direct, indirect, incidental, special,        #
# exemplary, or consequential damages (including, but not limited to,          #
# procurement of substitute goods or services; loss of use, data, or profits;  #
# or business interruption) however caused and on any theory of liability,     #
# whether in contract, strict liability, or tort (including negligence or      #
# otherwise) arising in any way out of the use of this software, even if       #
# advised of the possibility of such damage.                                   #
#                                                                              #
# You should have received a copy of the GNU General Public License            #
# along with IFC. If not, see <http://www.gnu.org/licenses/>.                  #
################################################################################

#' @title IFC_pops Tree
#' @description
#' Helper to extract hierarchy tree from pops.
#' @param pops list of tree characters.
#' @keywords internal
popsTree <- function(pops) {
  assert(pops, cla = c("IFC_pops","Affiliated"))
  
  # temp fun to reorder pops
  f <- function(P) {
    i=1;
    v = sapply(P, FUN = function(p) switch(p$type, "G"=1,"B"=2,"T"=3,"C"=4))
    P=P[unlist(use.names = F,by(cbind(names(v), v), v, FUN = function(d) sort(d[,1], decreasing = identical(head(d[,2],1),"1"))))]
    while(i<length(P)) {
      if(identical(P[[i]]$type,"G")) {
        index=which(P[[i]]$base==names(P))
        index=index[index>i]
        if(length(index)!=0) {
          P = c(P[setdiff(seq_len(index), i)],P[i],tail(P,-index))
        } else {
          i=i+1
        }
      } else {
        i=i+1
      }
    }
    class(P) = c("IFC_pops","Affiliated","Ordered")
    return(P)
  }
  
  # reorder pops having "B" 1st, "T" 2nd, "C" 3rd, ordered by name for each group
  # then "G" are positioned as sub elements of the former 3 also ordered by name
  P = f(pops)
  # extract level of each pop
  lev = popsGetLevels(P, min)
  
  # compute indentation (tree)
  indent = lapply(seq_along(lev), FUN = function(i) {
    p = P[[i]]
    if(p$name=="All") return("")
    v = c(rep(c("\u2502","\uffa0"),max(0,lev[i]-2)))#"│ﾠ"
    if(i == length(lev)) return(c(v,"\u2514"))      #"└"
    if(lev[i]  > lev[i + 1]) return(c(v,"\u2514"))  #"└"
    if(lev[i] == lev[i + 1]) return(c(v,"\u251c"))  #"├"
    if(lev[i]  < lev[i + 1]) return(c(v,"\u255e"))  #"╞"
  })
  
  # extract sub nodes indices
  sub_nodes = which(sapply(names(lev), FUN = function(i) P[[i]]$type != "G"))
  k = cbind(beg=sub_nodes, end=c(sub_nodes[-1], 1+length(lev)))
  
  # create matrix of characters
  m = max(sapply(indent, length))
  D = do.call(rbind, sapply(indent, simplify = FALSE, FUN = function(x) c(x, rep("", m - length(x)))))
  
  # finalize tree replacing final nodes "├" by "└"  or "╞" by "╘" and clearing empty branches "|"
  structure(as.list(unlist(use.names = FALSE, recursive = FALSE, lapply(seq_len(nrow(k)), FUN = function(i) {
    l = diff(k[i,]); S = seq_len(l); K = S + k[i,1] - 1;
    d = D[K, , drop = FALSE]
    n = (lev[rownames(k)[i]] - 1) * 2 + 1
    if(i == nrow(k)) n = 1
    while(n) {
      tmp = max(which(d[, n] != "\u2502"))
      if(tmp != l) {
        d[tmp, n] <- switch(d[tmp, n], "\u251c"="\u2514", "\u255e"="\u2558",d[tmp, n])
        d[tail(S, -tmp), n] <- "\uffa0"
        n = n + 2
      } else {
        n = 0
      }
    }
    apply(d, 1, paste0, collapse="")
  }))), names = names(lev))
}
