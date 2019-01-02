%{
/*============================================================================
 * Define the grammar for the mathematical expression
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2019 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_printf.h"

#include "mei_node.h"
#include "mei_parser_glob.h"

%}

/*----------------------------------------------------------------------------
 * %pure_parser ?
 *----------------------------------------------------------------------------*/

%union {
    double iValue;              /* double value */
    char sIndex[200];           /* variable, constant or function identifier */
    mei_node_t *nPtr;           /* node pointer */
};

%token <iValue> NUMBER
%token <sIndex> VAR FUN1 FUN2 FUN3 FUN4 IN1D
%token WHILE IF PRINT
%nonassoc IFX
%nonassoc ELSE

%right '='
%left OR
%left AND
%left EQ NE
%left GE LE '>' '<'
%left '+' '-'
%left '*' '/'
%right '!'
%nonassoc UMINUS
%nonassoc UPLUS
%right '^'

%type <nPtr> stmt expr stmt_list


%start program

%%

/*----------------------------------------------------------------------------
 * Yacc grammars
 *----------------------------------------------------------------------------*/

program: /* empty */                                   { return 1; }
        | stmt_list                                    { mei_glob_root = $1; return 1; }
        ;

stmt_list:
          stmt_list stmt                               { $$ = mei_opr_node(';', 2, $1, $2); }
        | stmt                                         { $$ = $1; }
        ;

stmt:
          ';'                                          { $$ = mei_opr_node(';', 2, NULL, NULL); }
        | expr ';'                                     { $$ = $1; }
        | PRINT expr ';'                               { $$ = mei_opr_node(PRINT, 1, $2); }
        | VAR '=' expr ';'                             { $$ = mei_opr_node('=', 2, mei_id_node($1), $3); }
        | WHILE '(' expr ')' stmt                      { $$ = mei_opr_node(WHILE, 2, $3, $5); }
        | IF '(' expr ')' stmt %prec IFX               { $$ = mei_opr_node(IF, 2, $3, $5); }
        | IF '(' expr ')' stmt ELSE stmt               { $$ = mei_opr_node(IF, 3, $3, $5, $7); }
        | '{' stmt_list '}'                            { $$ = $2; }
        | stmt error                                   { yyerror(mei_label_node($1)); return 0; }
        ;

expr:
          NUMBER                                       { $$ = mei_const_node($1); }
        | VAR                                          { $$ = mei_id_node($1); }
        | FUN1 '(' expr ')'                            { $$ = mei_func_node($1, $3); }
        | FUN2 '(' expr ',' expr ')'                   { $$ = mei_funcx_node($1, 2, $3, $5); }
        | FUN3 '(' expr ',' expr ',' expr ')'          { $$ = mei_funcx_node($1, 3, $3, $5, $7); }
        | FUN4 '(' expr ',' expr ',' expr ',' expr ')' { $$ = mei_funcx_node($1, 4, $3, $5, $7, $9); }
        | '!' expr                                     { $$ = mei_opr_node('!',     1, $2); }
        | '-' expr %prec UMINUS                        { $$ = mei_opr_node(UMINUS,  1, $2); }
        | '+' expr %prec UPLUS                         { $$ = mei_opr_node(UPLUS,   1, $2); }
        | expr '+' expr                                { $$ = mei_opr_node('+', 2, $1, $3); }
        | expr '-' expr                                { $$ = mei_opr_node('-', 2, $1, $3); }
        | expr '*' expr                                { $$ = mei_opr_node('*', 2, $1, $3); }
        | expr '/' expr                                { $$ = mei_opr_node('/', 2, $1, $3); }
        | expr '^' expr                                { $$ = mei_opr_node('^', 2, $1, $3); }
        | expr '<' expr                                { $$ = mei_opr_node('<', 2, $1, $3); }
        | expr '>' expr                                { $$ = mei_opr_node('>', 2, $1, $3); }
        | expr GE expr                                 { $$ = mei_opr_node(GE,  2, $1, $3); }
        | expr LE expr                                 { $$ = mei_opr_node(LE,  2, $1, $3); }
        | expr NE expr                                 { $$ = mei_opr_node(NE,  2, $1, $3); }
        | expr EQ expr                                 { $$ = mei_opr_node(EQ,  2, $1, $3); }
        | expr OR expr                                 { $$ = mei_opr_node(OR,  2, $1, $3); }
        | expr AND expr                                { $$ = mei_opr_node(AND, 2, $1, $3); }
        | '(' expr ')'                                 { $$ = $2; }
        | expr error                                   { yyerror(mei_label_node($1)); return 0; }
        ;

%%

/*----------------------------------------------------------------------------
 * Parsing error management
 *----------------------------------------------------------------------------*/

void
yyerror(const char *s)
{
    int l;

    mei_free_node(mei_glob_root);
    --mei_glob_column;

    /* bft_printf("Warning: %s\n", s); */
    /* bft_printf("line %i column %i \n", line, column); */

    BFT_REALLOC(mei_glob_label_list,  mei_glob_ierr_list+1, char*);
    BFT_REALLOC(mei_glob_line_list,   mei_glob_ierr_list+1, int);
    BFT_REALLOC(mei_glob_column_list, mei_glob_ierr_list+1, int);

    l = strlen("Warning: ") +1;
    BFT_MALLOC(mei_glob_label_list[mei_glob_ierr_list], l, char);
    strncpy(mei_glob_label_list[mei_glob_ierr_list], "Error: ", l);

    l += strlen(s);
    BFT_REALLOC(mei_glob_label_list[mei_glob_ierr_list], l, char);
    strncat(mei_glob_label_list[mei_glob_ierr_list], s, l);

    l += strlen(" \n");
    BFT_REALLOC(mei_glob_label_list[mei_glob_ierr_list], l, char);
    strncat(mei_glob_label_list[mei_glob_ierr_list], " \n", l);

    mei_glob_line_list[mei_glob_ierr_list]   = mei_glob_line;
    mei_glob_column_list[mei_glob_ierr_list] = mei_glob_column;

    mei_glob_ierr_list++;

    return;
}

/*----------------------------------------------------------------------------*/
