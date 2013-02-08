#include <Python.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

typedef struct {
  double a;
  double c;
  double g;
  double t;
} column;

typedef struct {
  int pos;
  double score;
} pbs;

int compare_pbs(const void *a, const void *b) {
  /* comparison function to sort putative binding sites */
  const pbs *pbsa = (const pbs *) a;
  const pbs *pbsb = (const pbs *) b;
  // by pos
  //return (pbsa->pos < pbsb->pos);
  // by score
  return (pbsa->score < pbsb->score);
}

column* create_PSSM(char **sites, size_t num_sites, column bg_prob) {
  /* Given list of sites and background probabilities, return the PSSM as a list of
  columns, each containing score for a,c,t,g.
  */
  int i, j;
  int site_length = strlen(sites[0]);
  column *PSSM = (column *) malloc(site_length * sizeof(column));
  column base_counts; // base counts for a column
  //printf("site length %d\n", site_length);
  //printf("numsites %d\n", num_sites);
  
  // count bases on each column
  for (i = 0; i < site_length; i++) {
    //base_counts.a = base_counts.c = base_counts.g = base_counts.t = 1e-10;
    // set pseudocounts
    base_counts.a = bg_prob.a;
    base_counts.c = bg_prob.c;
    base_counts.g = bg_prob.g;
    base_counts.t = bg_prob.t;
    
    //printf("%d %d %d %d\n", base_counts.a, base_counts.c, base_counts.g, base_counts.t);
    for (j = 0; j < num_sites; ++j) {
      switch (sites[j][i]) {
      case 'A':
      case 'a': base_counts.a += 1; break;
      case 'C':
      case 'c': base_counts.c += 1; break;
      case 'G':
      case 'g': base_counts.g += 1; break;
      case 'T':
      case 't': base_counts.t += 1; break;
      }
    }
    // normalize counts
    int total = (base_counts.a + base_counts.c + base_counts.g + base_counts.t);
    base_counts.a /= total;
    base_counts.c /= total;
    base_counts.g /= total;
    base_counts.t /= total;
    // fill PSSM column
    PSSM[i].a = log2(base_counts.a / bg_prob.a);
    PSSM[i].c = log2(base_counts.c / bg_prob.c);
    PSSM[i].g = log2(base_counts.g / bg_prob.g);
    PSSM[i].t = log2(base_counts.t / bg_prob.t);
  }
  return PSSM;
}

double score_seq(column *PSSM, int PSSM_length, char *seq) {
  /* Given PSSM and a sequence (of same length with PSSM, return score */
  int i;
  double score = 0.0;
  for (i = 0; i < PSSM_length; ++i) {
    switch (seq[i]) {
    case 'A':
    case 'a': score += PSSM[i].a; break;
    case 'C':
    case 'c': score += PSSM[i].c; break;
    case 'G':
    case 'g': score += PSSM[i].g; break;
    case 'T':
    case 't': score += PSSM[i].t; break;
    }
  }
  return score;
}

pbs* site_search(char **sites, int num_sites, char* seq, column bg_prob, int *num_pbs) {
  /* Build PSSM using sites, and scan the genome using the PSSM and return the
     list of putative binding sites sorted by score.
  */
  int i;
  double score;
  int num_possible_pbs = strlen(seq) - strlen(sites[0]) + 1;
  pbs *putative_binding_sites = (pbs *) malloc(num_possible_pbs * sizeof(pbs));
  column* PSSM = create_PSSM(sites, num_sites, bg_prob);

  *num_pbs = 0;
  for (i = 0; i < num_possible_pbs; i++) {
    score = score_seq(PSSM, strlen(sites[0]), seq+i);
    if (score > 0.0) {
      putative_binding_sites[*num_pbs].pos = i;
      putative_binding_sites[*num_pbs].score = score;
      *num_pbs += 1;
    }
  }
  free(PSSM);
  return putative_binding_sites;
}

/*
 * Function to be called from Python
 */
static PyObject* yassi_search(PyObject* self, PyObject* args) {
  char *seq;
  char **sites;
  int num_sites;
  pbs *putative_binding_sites;
  int num_putative_binding_sites;
  int i;
  PyObject *list_obj; // to parse list of sites (strings)
  PyObject *bg_prob_obj = NULL;
  PyObject *return_list, *tuple;
  column bg_prob;
  
  if (!PyArg_ParseTuple(args, "O!s|O!", &PyList_Type, &list_obj, &seq, &PyList_Type, &bg_prob_obj))
    return NULL;
  
  /* set background probabilities */
  if (bg_prob_obj) {
    bg_prob.a = PyFloat_AsDouble(PyList_GetItem(bg_prob_obj, 0));
    bg_prob.c = PyFloat_AsDouble(PyList_GetItem(bg_prob_obj, 1));
    bg_prob.g = PyFloat_AsDouble(PyList_GetItem(bg_prob_obj, 2));
    bg_prob.t = PyFloat_AsDouble(PyList_GetItem(bg_prob_obj, 3));
  } else {
    bg_prob.a = bg_prob.c = bg_prob.g = bg_prob.t = 0.25;
  }

  // get number of strings passed
  num_sites = PyList_Size(list_obj);
  sites = (char **) malloc(num_sites * sizeof(char*));
  for (i = 0; i < num_sites; ++i) {
    sites[i] = PyString_AsString(PyList_GetItem(list_obj, i));
  }
  // find sites
  putative_binding_sites = site_search(sites, num_sites, seq, bg_prob, &num_putative_binding_sites);
  // sort putative binding sites by scores
  qsort(putative_binding_sites, num_putative_binding_sites, sizeof(pbs), compare_pbs);
  
  // pass the list of (pos, score) back to Python
  return_list = PyList_New(num_putative_binding_sites);
  for (i = 0; i < num_putative_binding_sites; ++i) {
    tuple = Py_BuildValue("id", putative_binding_sites[i].pos, putative_binding_sites[i].score);
    PyList_SetItem(return_list, i, tuple);
  }
  free(putative_binding_sites);
  free(sites);
  PyObject *result = Py_BuildValue("O", return_list);
  Py_DECREF(return_list);
  return result;
}

static PyObject* yassi_build_PSSM(PyObject* self, PyObject* args) {
  PyObject *list_obj, *bg_prob_obj = NULL;
  int num_sites;
  char **sites;
  int i;
  column bg_prob;
  column *PSSM;
  int PSSM_len;
  PyObject *col;
  PyObject *return_PSSM;
  
  if (!PyArg_ParseTuple(args, "O!|O!", &PyList_Type, &list_obj, &PyList_Type, &bg_prob_obj))
    return NULL;

  /* set background probabilities */
  if (bg_prob_obj) {
    bg_prob.a = PyFloat_AsDouble(PyList_GetItem(bg_prob_obj, 0));
    bg_prob.c = PyFloat_AsDouble(PyList_GetItem(bg_prob_obj, 1));
    bg_prob.g = PyFloat_AsDouble(PyList_GetItem(bg_prob_obj, 2));
    bg_prob.t = PyFloat_AsDouble(PyList_GetItem(bg_prob_obj, 3));
  } else {
    bg_prob.a = bg_prob.c = bg_prob.g = bg_prob.t = 0.25;
  }
  
  num_sites = PyList_Size(list_obj);
  sites = (char **) malloc(num_sites * sizeof(char*));
  for (i = 0; i < num_sites; ++i) {
    sites[i] = PyString_AsString(PyList_GetItem(list_obj, i));
  }
  
  PSSM = create_PSSM(sites, num_sites, bg_prob);
  PSSM_len = strlen(sites[0]);
  // return PSSM as Python object (list of ACTG weights, in order)
  return_PSSM = PyList_New(PSSM_len);
  for (i = 0; i < PSSM_len; ++i) {
    col = Py_BuildValue("dddd", PSSM[i].a, PSSM[i].c, PSSM[i].g, PSSM[i].t);
    PyList_SetItem(return_PSSM, i, col);
  }
  return Py_BuildValue("O", return_PSSM);
}

/*
 * Bind Python function names to our C functions
 */
static PyMethodDef yassi_methods[] = {
  {"search", yassi_search, METH_VARARGS},
  {"build_PSSM", yassi_build_PSSM, METH_VARARGS},
  {NULL, NULL, 0}
};

/*
 * Python calls this to let us initialize our module
 */
void inityassi() {
  (void) Py_InitModule("yassi", yassi_methods);
}
