/* source code courtesy of Frank Thomas Braun */

/* calc3d.c: Generation of the graph of the syntax tree */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mei_node.h"
#include "mei_parser.h"

int del = 1; /* distance of graph columns */
int eps = 3; /* distance of graph lines */

/* interface for drawing (can be replaced by "real" graphic using GD or other) */
int  graphik(mei_node_t *p);
void graphInit(void);
void graphFinish(void);
void graphBox(const char *s, int *w, int *h);
void graphDrawBox(const char *s, int c, int l);
void graphDrawArrow(int c1, int l1, int c2, int l2);
void graphTest(int l, int c);

/* recursive drawing of the syntax tree */
void exNode(mei_node_t *p, int c, int l, int *ce, int *cm);

/*****************************************************************************/

/* main entry point of the manipulation of the syntax tree */
int
graphik(mei_node_t *p)
{
    int rte, rtm;

    graphInit();
    exNode(p, 0, 0, &rte, &rtm);
    graphFinish();
    return 0;
}

/*c----cm---ce---->                       drawing of leaf-nodes
 l leaf-info
 */

/*c---------------cm--------------ce----> drawing of non-leaf-nodes
 l            node-info
 *                |
 *    -------------     ...----
 *    |       |               |
 *    v       v               v
 * child1  child2  ...     child-n
 *        che     che             che
 *cs      cs      cs              cs
 *
 */

void
exNode(mei_node_t *p,
       int c, int l,        /* start column and line of node */
       int *ce, int *cm)    /* resulting end column and mid of node */
{
    int w, h;         /* node width and height */
    const char *s;    /* node text */
    int cbar;         /* "real" start column of node (centred above subnodes) */
    int k;            /* child number */
    int che, chm;     /* end column and mid of children */
    int cs;           /* start column of children */
    char word[20];    /* extended node text */

    if (!p) {
      *ce = c;
      *cm = c;
      return;
    }

    strcpy(word, "???"); /* should never appear */
    s = word;
    switch(p->flag) {
        case FUNC1: strcpy(word, p->type->func.name); break;
        case FUNC2: strcpy(word, p->type->funcx.name); break;
        case FUNC3: break;
        case FUNC4: break;
        case CONSTANT: sprintf(word, "c(%f)", p->type->con.value); break;
        case ID:  sprintf(word, "id(%s)", p->type->id.i); break;
        case OPR:
            switch(p->type->opr.oper){
                case WHILE:     s = "while"; break;
                case IF:        s = "if";    break;
                case PRINT:     s = "print"; break;
                case ';':       s = "[;]";     break;
                case '=':       s = "[=]";     break;
                case UMINUS:    s = "[-]";     break;
                case UPLUS :    s = "[+]";     break;
                case '+':       s = "[+]";     break;
                case '-':       s = "[-]";     break;
                case '*':       s = "[*]";     break;
                case '/':       s = "[/]";     break;
                case '^':       s = "[^]";     break;
                case '<':       s = "[<]";     break;
                case '>':       s = "[>]";     break;
                case '!':       s = "[!]";     break;
                case GE:        s = "[>=]";    break;
                case LE:        s = "[<=]";    break;
                case NE:        s = "[!=]";    break;
                case EQ:        s = "[==]";    break;
                case OR:        s = "[||]";    break;
                case AND:       s = "[&&]";    break;
            }
            break;
    }

    /* construct node text box */
    graphBox(s, &w, &h);
    cbar = c;
    *ce = c + w;
    *cm = c + w / 2;

    /* node is leaf */
    if (p->flag == CONSTANT || p->flag == ID) {
        graphDrawBox(s, cbar, l);
        return;
    }

    /* node has children */
    cs = c;
    if (p->flag == FUNC1) {
      exNode(p->type->func.op, cs, l+h+eps, &che, &chm);
      cs = che;
    } else if (p->flag == FUNC2) {
      for (k = 0; k < p->type->funcx.nops; k++) {
        exNode(p->type->funcx.op[k], cs, l+h+eps, &che, &chm);
        cs = che;
      }

    } else {
      for (k = 0; k < p->type->opr.nops; k++) {
          exNode(p->type->opr.op[k], cs, l+h+eps, &che, &chm);
          cs = che;
      }
    }

    /* total node width */
    if (w < che - c) {
        cbar += (che - c - w) / 2;
        *ce = che;
        *cm = (c + che) / 2;
    }

    /* draw node */
    graphDrawBox(s, cbar, l);

    /* draw arrows (not optimal: children are drawn a second time) */
    cs = c;
    if (p->flag == FUNC1) {
      exNode(p->type->func.op, cs, l+h+eps, &che, &chm);
      graphDrawArrow(*cm, l+h, chm, l+h+eps-1);
      cs = che;
    } else if (p->flag == FUNC2) {
      for (k = 0; k < p->type->funcx.nops; k++) {
          exNode(p->type->funcx.op[k], cs, l+h+eps, &che, &chm);
          graphDrawArrow(*cm, l+h, chm, l+h+eps-1);
          cs = che;
      }
    } else {
      for (k = 0; k < p->type->opr.nops; k++) {
          exNode(p->type->opr.op[k], cs, l+h+eps, &che, &chm);
          graphDrawArrow(*cm, l+h, chm, l+h+eps-1);
          cs = che;
      }
    }
}

/* interface for drawing */

#define lmax 200
#define cmax 200

char graph[lmax][cmax]; /* array for ASCII-Graphic */
int graphNumber = 0;

void
graphTest(int l,
          int c)
{   int ok;
    ok = 1;
    if (l < 0) ok = 0;
    if (l >= lmax) ok = 0;
    if (c < 0) ok = 0;
    if (c >= cmax) ok = 0;
    if (ok) return;
    printf("\n+++error: l=%d, c=%d not in drawing rectangle 0, 0 ... %d, %d",
        l, c, lmax, cmax);
    exit(1);
}

void
graphInit(void)
{
    int i, j;
    for (i = 0; i < lmax; i++) {
        for (j = 0; j < cmax; j++) {
            graph[i][j] = ' ';
        }
    }
}

void
graphFinish(void)
{
    int i, j;
    for (i = 0; i < lmax; i++) {
        for (j = cmax-1; j > 0 && graph[i][j] == ' '; j--);
        graph[i][cmax-1] = 0;
        if (j < cmax-1) graph[i][j+1] = 0;
        if (graph[i][j] == ' ') graph[i][j] = 0;
    }
    for (i = lmax-1; i > 0 && graph[i][0] == 0; i--);
    printf("\n\nGraph %d:\n", graphNumber++);
    for (j = 0; j <= i; j++) printf("\n%s", graph[j]);
    printf("\n");
}

void
graphBox(const char *s,
         int        *w,
         int        *h)
{
    *w = strlen(s) + del;
    *h = 1;
}

void
graphDrawBox(const char *s,
             int         c,
             int         l)
{
    int i;
    graphTest(l, c+strlen(s)-1+del);
    for (i = 0; (unsigned) i < strlen(s); i++) {
        graph[l][c+i+del] = s[i];
    }
}

void
graphDrawArrow(int c1,
               int l1,
               int c2,
               int l2)
{
    int m;
    graphTest(l1, c1);
    graphTest(l2, c2);
    m = (l1 + l2) / 2;
    while (l1 != m) { graph[l1][c1] = '|'; if (l1 < l2) l1++; else l1--; }
    while (c1 != c2) { graph[l1][c1] = '-'; if (c1 < c2) c1++; else c1--; }
    while (l1 != l2) { graph[l1][c1] = '|'; if (l1 < l2) l1++; else l1--; }
    graph[l1][c1] = '|';
}

