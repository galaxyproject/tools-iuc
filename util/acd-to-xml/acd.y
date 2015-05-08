%{
#include <stdio.h>
#include <string.h>
extern "C" int yyparse (void);
extern int yylex (void);

void yyerror(const char *str)
{
    fprintf(stderr,"error: %s\n",str);
}

int yywrap()
{
    return 1;
}

int main()
{
    yyparse();
}

%}

%union {
  char* str;
}

%token APPTOK QUOTE OBRACE EBRACE OBRACKET EBRACKET COLON OSEC ESEC;
%token <str> QWORDS WORDS WORD;

%%
application:
    application_section
    sections
    {
        printf("</application>\n");
    }
;

application_section:
    APPTOK COLON WORD 
    {
        printf("<?xml version=\"1.0\"?>\n<application id=\"%s\">\n", $3);
    }
    OBRACKET
    { printf("\t<metadata>\n"); }
    metadata_list
    { printf("\t</metadata>\n"); }
    EBRACKET
;

sections: %empty
    | sections section;

section:
    OSEC COLON WORD
    {
        printf("\t<section id=\"%s\">\n", $3);
    }
    OBRACKET
    { printf("\t<metadata>\n"); }
    metadata_list
    { printf("\t</metadata>\n"); }
    EBRACKET
    parameters_and_sections
    ESEC COLON WORD
    {
        printf("\t</section>\n");
    }
;

parameters_and_sections: %empty
    | parameters_and_sections parameter_list
    | parameters_and_sections section ;

parameter_list: %empty
              | parameter_list parameter;

parameter:
    WORD COLON WORD
    { printf("\t<parameter type=\"%s\" name=\"%s\">\n", $1, $3);}
    OBRACKET
    metadata_list
    EBRACKET
    { printf("</parameter>");}

metadata_list: %empty
        | metadata_list metadata;

metadata: WORD COLON QWORDS
        {
            printf("\t\t<%s><![CDATA[%s]]></%s>\n", $1, $3, $1);
        }
        ;
%%
