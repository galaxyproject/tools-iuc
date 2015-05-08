%{
#include <stdio.h>
#include <string.h>
#include "y.tab.h"
%}

%%
#.*\n                   /* ignore comments */;
application             return APPTOK;
section                 return OSEC;
endsection              return ESEC;
variable                return VARIABLE;
\"[^\"]*\"              yylval.str=strndup(yytext + 1, strlen(yytext) - 2); return QWORDS;
[a-zA-Z][a-zA-Z0-9]*    yylval.str=strdup(yytext); return WORD;
\"                      return QUOTE;
\{                      return OBRACE;
\}                      return EBRACE;
\[                      return OBRACKET;
\]                      return EBRACKET;
:                       return COLON;
\n                      /* ignore EOL */;
[ \t]+                  /* ignore whitespace */;
%%
