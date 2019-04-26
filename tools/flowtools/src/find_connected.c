#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

//static const char *rcsid = "$Id: find_connected.c,v 1.1 2008/09/05 21:54:40 rpl Exp $";

int find_connected(int **G, int num_dense_grids, int ndim, int *grid_clusterID);
void depth_first(int startnode);

void bail(const char *);        /* exits via abort */
static void check_clusters(int *gcID, int ndense);
static void merge_cluster(int from, int into);
  

/* Vars that will not change througout the depth-first recursion.  We
   store them here to avoid endless replication on the stack. */
static int **Gr=0;
static int *gcID = 0;          /* grid cluster IDs */
static int *cluster_count=0;   /* count of nodes per cluster */
static int ndense=0;
static int ndim=0;
/* cid changes between depth-first searches, but is constant within a
   single search, so it goes here. */
static int cid=0;

/* Find connected components in the graph of neighboring grids defined in G.
 *
 * Output:
 *
 * grid_clusterID[]  -- cluster to which each dense grid was assigned
 * return value      -- number of clusters assigned.
 */
int find_connected(int **G, int n_dense_grids, int num_dm, int *grid_clusterID)
{
  int nclust=0;                  /* number of clusters found */
  int i;
  int *subfac;
  int subval=0,nempty=0;
    int clustid=0;
  
  size_t sz = n_dense_grids*sizeof(int); 
  cluster_count = malloc(sz);
  if(!cluster_count)
    bail("find_connected: Unable to allocate %zd bytes.\n");
  memset(cluster_count,0,sz);

  /* set up the statics that will be used in the DFS */
  Gr=G;
  gcID = grid_clusterID;
  ndense = n_dense_grids;
  ndim = num_dm;

  
  for(i=0;i<ndense;++i)
    grid_clusterID[i] = -1;

  for(i=0;i<ndense;++i) {
    if(grid_clusterID[i] < 0) {  /* grid hasn't been assigned yet */
      cid = nclust++;
      depth_first(i);
    }
  }

#ifndef NDEBUG
  check_clusters(gcID,ndense);
#endif 
  
  /* At this point we probably have some clusters that are empty due to merging.
     We want to compact the cluster numbering to eliminate the empty clusters. */

  subfac = malloc(sz);
  if(!subfac)
    bail("find_connected: Unable to allocate %zd bytes.\n");
 subval=0;
 nempty=0;

  /* cluster #i needs to have its ID decremented by 1 for each empty cluster with
     ID < i.  Precaclulate the decrements in this loop: */
  for(i=0;i<nclust;++i) {
    //clustid = grid_clusterID[i];
    if(cluster_count[i] == 0) {
      subval++;
      nempty++;
    }
    subfac[i] = subval;
  }

  /* Now apply the decrements to all of the dense grids */
  for(i=0;i<ndense;++i) {
   clustid = grid_clusterID[i];
    grid_clusterID[i] -= subfac[clustid];
  }

#ifndef NDEBUG
  //  check_clusters(gcID,ndense);
#endif  
  
  /* correct the number of clusters found */
  nclust -= nempty;

  return nclust;
}


/* Do a depth-first search for a single connected component in graph
 * G.  Start from node, tag the nodes found with cid, and record
 * the tags in grid_clusterID.  Also, record the node count in
 * cluster_count.  If we find a node that has already been assigned to
 * a cluster, that means we're merging two clusters, so zero out the
 * old cid's node count.
 *
 * Note that our graph is constructed as a DAG, so we can be
 * guaranteed to terminate without checking for cycles.  
 *
 * Note2: this function can potentially recurse to depth = ndense.
 * Plan stack size accordingly.  
 *
 * Output:
 *
 * grid_clusterID[] -- array where we tag the nodes.
 * cluster_count[]  -- count of the number of nodes per cluster.
 */
   
void depth_first(int node)
{
  int i;

  if(gcID[node] == cid)         // we're guaranteed no cycles, but it is possible to reach a node 
    return;                     // through two different paths in the same cluster.  This early
                                // return saves us some unneeded work.

  /* Check to see if we are merging a cluster */
  if(gcID[node] >= 0) {
    /* We are, so zero the count for the old cluster. */
    cluster_count[ gcID[node] ] = 0;  
    merge_cluster(gcID[node], cid);
    return;
  }

  /* Update for this node */
  gcID[node] = cid;
  cluster_count[cid]++;

  /* Recursively search the child nodes */
  for(i=0; i<ndim; ++i)
    if(Gr[node][i] >= 0)      /* This is a child node */
      depth_first(Gr[node][i]);
}

void bail(const char *msg)
{
  fprintf(stderr,"%s",msg);
  abort();
}


static void check_clusters(int *gcID, int ndense)
{
  int i;

  for(i=0; i<ndense; ++i)
    if(gcID[i] < 0) {
      fprintf(stderr,"faulty cluster id at i= %d\n",i);
      abort();
    }
}

static void merge_cluster(int from, int into)
{
  int i;

  for(i=0; i<ndense; ++i)
    if(gcID[i] == from)
      gcID[i] = into;
}
