%{
/*============================================================================
 * Define the scanner for the mathematical expression
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "bft_mem.h"
#include "bft_error.h"

#include "mei_hash_table.h"
#include "mei_evaluate.h"
#include "mei_node.h"
#include "mei_parser.h"
#include "mei_parser_glob.h"

/* Redefine macro to redirect scanner input from string instead of standard input. */

#undef YY_INPUT
#define YY_INPUT(b, r, ms) ( r = my_yyinput(b, ms) )

int my_yyinput(char *buffer, int max_size);

char *buff = NULL;
size_t len;

%}

%%
" "         {mei_glob_column++;                  } /* skip whitespace */
"\t"        {mei_glob_column+=8;                 } /* skip whitespace */
"\n"        {mei_glob_line++; mei_glob_column=0; } /* skip whitespace */
"#".*\n     {mei_glob_line++; mei_glob_column=0; } /* discard commentary (other possibility: "#"[^\n]*"\n" ) */


([0-9]+|[0-9]+\.[0-9]*|[0-9]*\.[0-9]+)([Ee][+\-]?[0-9]+)? {
    yylval.iValue = strtod(yytext, NULL);
    mei_glob_column += yyleng;
    return NUMBER;
}


"e"|"pi" {
    strncpy(yylval.sIndex, yytext, yyleng+1);
    mei_glob_column += yyleng;
    return VAR;
}

"min"|"max"|"atan2"|"mod" {
    strncpy(yylval.sIndex, yytext, yyleng+1);
    mei_glob_column += yyleng;
    return FUN2;
}

"exp"|"log"|"sqrt"|"sin"|"cos"|"tan"|"asin"|"acos"|"atan"|"sinh"|"cosh"|"tanh"|"abs"|"int" {
    strncpy(yylval.sIndex, yytext, yyleng+1);
    mei_glob_column += yyleng;
    return FUN1;
}

"while"         {mei_glob_column+=5; return WHILE; }
"if"            {mei_glob_column+=2; return IF;    }
"else"          {mei_glob_column+=4; return ELSE;  }
"print"         {mei_glob_column+=5; return PRINT; }


[_a-zA-Z][\[\]\._a-zA-Z0-9]* {
    strncpy(yylval.sIndex, yytext, yyleng+1);
    if (yyleng+1 > 200)
        bft_error(__FILE__, __LINE__, 0, "Error: identifier is to long\n");
    mei_glob_column += yyleng;
    return VAR;
}


">="             {mei_glob_column+=2; return GE;  }
"<="             {mei_glob_column+=2; return LE;  }
"=="             {mei_glob_column+=2; return EQ;  }
"!="             {mei_glob_column+=2; return NE;  }
"||"             {mei_glob_column+=2; return OR;  }
"&&"             {mei_glob_column+=2; return AND; }
"!"              {mei_glob_column++;  return *yytext;}


[-()<>=+^*/;{},]  {mei_glob_column++; return *yytext;}


. {
    mei_glob_column++;
    len = strlen("unknown character: ") +1;
    BFT_MALLOC(buff, len, char);
    strncpy(buff, "unknown character: ", len);
    BFT_REALLOC(buff, len+yyleng, char);
    strncat(buff, yytext, yyleng);
    yyerror(buff);
    BFT_FREE(buff);
}

%%

int yywrap(void)
{
    return 1;
}


/* The calling sequence for 'YY_INPUT()' is 'YY_INPUT(buff,result,max_size)'.
   Its action is to place up to 'max_size' characters in the character
   array 'buff' and return in the integer variable 'buffer_length' either
   the number of characters read or the constant 'YY_NULL' to indicate 'EOF'.
   The default 'YY_INPUT' reads from the global file-pointer 'yyin'. */

/* max_sise = 8192 bytes. One can change this default
   value with YY_BUF_SIZE macro */


int
my_yyinput(char *buffer, int max_size)
{
    int buffer_length;

    if (max_size < mei_glob_string_end - mei_glob_string_begin)
    {
        buffer_length = max_size;
    }
    else
    {
        buffer_length = mei_glob_string_end - mei_glob_string_begin;
    }

    if (buffer_length > 0)
    {
        memcpy(buffer, mei_glob_string_begin, buffer_length);
        mei_glob_string_begin +=buffer_length;
    }

    return buffer_length;
}
