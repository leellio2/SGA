/******************************************************************************
 Purpose:
    This macro can be used as diagnostic tool in propensity score weighting
    analyses with subgroups. The macro produce the following output:
    * Tables with descriptive statistics and unweighted and weighted standardized
    * differences of confounders by subgroup
    * connectS plots to display standardized differences
    * Histograms of propensity score/prognostic variables by treatment group

 Author: Daniel Wojdyla
 Date: August 2022
 Version: 0.9

 See documentation for description of parameters
 The stddiff macro is an adapted version from the macro available at:
    https://www.lerner.ccf.org/qhs/software/lib/stddiff.sas
*********************************************************************************/

* ######################################################################################### *
* ################################### diagSGA macro ####################################### *
* ######################################################################################### *;
%macro diagSGA(ds=,weights=,propens=,w_type=O,trtvar=,contconf=,catconf=,subgrps=,prognos=,
               outdir=,name=analysis,conffmt=,subgfmt=,trtfmt=,vi=Y,ss=Y,hgt=,wdt=,
               reflines=,fsize=,csize=,nrows=,ncols=,color=N,xmin=,xmax=,xstep=,hxlabel=);
   ods listing close;
   data _temp_;
      set &ds;

      %if %length(&weights.) > 0 %then %do;
         wvar = &weights;
         %end;

      %else %do;
         %if &w_type = I %then %do;
            wvar = (&trtvar=1)*(1/&propens.) + (&trtvar=0)*(1/(1-&propens.));
            %end;
         %else %do;
            wvar = (&trtvar=1)*(1-&propens.) + (&trtvar=0)*&propens.;
            %end;
         %end;
   run;

   * remove multiple blanks and get the total number of subgroups;
	%let subgvars = %sysfunc(compbl(&subgrps.));
	%let subgvars = %sysfunc(strip(&subgvars.));
	%let n_subgrps = %sysfunc(countc(&subgvars.,' '));
	%let n_subgrps = %eval(&n_subgrps. + 1);

	%do yy = 1 %to &n_subgrps.;
    	%let subgvar = %sysfunc(scan(&subgvars.,&yy.,' '));
         proc sort data=_temp_ out=sorted;
            by &subgvar.;
         run;

         data _temp2_;
            set sorted;
            by &subgvar;

            retain sgn 0;
            if first.&subgvar then do;
               sgn + 1;
               end;

            subg1 = &&yy;
            subg2 = sgn;
         run;

         proc append base=stacked data=_temp2_;
         run;

         proc sort data=_temp2_ out=_temp3_(keep=&subgvar) nodupkey;
            by &subgvar;
         run;

         data _temp3_;
            attrib var length=$32.;
            set _temp3_;
            var = "&subgvar";
            rename &subgvar=value;
         run;

         proc append base=labels data=_temp3_;
         run;
      %end;

   data labels;
      attrib label length=$25.;
      set labels end=eof;
      retain start;
      retain fmtname 'subgbf' type 'n';
      start + 1;
      label = compress(var) || "=" || compress(put(value,2.));
      output;

      if eof=1 then do;
         start + 1;
         label = "Overall";
         output;
         end;
         keep label start fmtname type;
      run;

   proc format library=work cntlin=labels;
   run;

   data stacked;
      set stacked
          _temp_(in=a);

      if a=1 then do;
         subg1 = &yy.;
         subg2 = 1;
         end;
   run;

   data stacked;
      set stacked;
      by subg1 subg2;

      retain subgn 0;
      if first.subg2 then subgn + 1;
   run;

   proc sql noprint;
	   select max(subgn) into: n_levels from stacked;
   quit;

   %do zz = 1 %to &n_levels.;
      data sd_data;
         set stacked(where=(subgn=&zz.));
      run;

      %stddiff(inds=sd_data,
               groupvar=&trtvar,
               numvars=&contconf,
               charvars=&catconf,
               wtvar=wvar,
               outds=results);

      * Standardized Differences for connectS plots;
      data for_plot;
         attrib standiff format=7.4;
         set results;

         subgn = &zz.;

         standiff = input(stddiff,7.4);
         if standiff ^= .;
         keep VarName subgn standiff;
      run;

      proc append base=stacked_sd data=for_plot;
      run;

      * Data for Standardized Differences Tables;
      data for_table;
         attrib standiff format=7.3;
         attrib standifft length=$10.;
         set results;

         subgn = &zz.;
         standiff = input(stddiff,7.3);
         if standiff ^= . then standifft=put(standiff,7.3);
         drop stddiff;
      run;

      proc means data=sd_data p5 nway;
         class &trtvar.;
         var &contconf;
         weight wvar;
         ods output Summary=extreme5;
      run;

      proc means data=sd_data p95 nway;
         class &trtvar.;
         var &contconf;
         weight wvar;
         ods output Summary=extreme95;
      run;

      proc means data=sd_data p50 nway;
         class &trtvar.;
         var &contconf;
         weight wvar;
         ods output Summary=extreme50;
      run;

      proc transpose data=extreme5(drop=NObs &trtvar.) out=extreme5_t;
      run;

      proc transpose data=extreme95(drop=NObs &trtvar.) out=extreme95_t;
      run;

      proc transpose data=extreme50(drop=NObs &trtvar.) out=extreme50_t;
      run;

      data percentiles;
         attrib perctls_w0 perctls_w1 length=$30.;
         merge extreme5_t(rename=(col1=p05_0 col2=p05_1))
               extreme95_t(rename=(col1=p95_0 col2=p95_1))
               extreme50_t(rename=(col1=p50_0 col2=p50_1));

         if p50_0 > 10 then perctls_w0 = compress(put(p50_0,7.)) || " [" || compress(put(p05_0,7.)) || " - " || compress(put(p95_0,7.)) || "]";
         else perctls_w0 = compress(put(p50_0,7.1)) || " [" || compress(put(p05_0,7.1)) || " - " || compress(put(p95_0,7.1)) || "]";

         if p50_1 > 10 then perctls_w1 = compress(put(p50_1,7.)) || " [" || compress(put(p05_1,7.)) || " - " || compress(put(p95_1,7.)) || "]";
         else perctls_w1 = compress(put(p50_1,7.1)) || " [" || compress(put(p05_1,7.1)) || " - " || compress(put(p95_1,7.1)) || "]";
         keep perctls_w0 perctls_w1;
      run;

      data for_table;
         merge for_table
               percentiles;
      run;

      proc append base=stacked_tab data=for_table;
      run;

      * Unweighted Versions;
      %stddiff(inds=sd_data,
               groupvar=&trtvar,
               numvars=&contconf,
               charvars=&catconf,
               outds=results_unw);


      * Standardized Differences for connectS plots (unweighted version);
      data for_plot_unw;
         attrib standiff_unw format=7.4;
         set results_unw;

         subgn = &zz.;

         standiff_unw = input(stddiff,7.4);
         if standiff_unw ^= .;
         keep VarName subgn standiff_unw;
      run;

      proc append base=stacked_sdunw data=for_plot_unw;
      run;

      * Data for Standardized Differences Tables (unweighted version);
      data for_table_unw;
         attrib standiff format=7.3;
         attrib standifft length=$10.;
         set results_unw;

         subgn = &zz.;
         standiff = input(stddiff,7.3);
         if standiff ^= . then standifft=put(standiff,7.3);
         drop stddiff;
      run;

      proc means data=sd_data p5 nway;
         class &trtvar;
         var &contconf;
         ods output Summary=extreme5U;
      run;

      proc means data=sd_data p95 nway;
         class &trtvar;
         var &contconf;
         ods output Summary=extreme95u;
      run;

      proc means data=sd_data p50 nway;
         class &trtvar;
         var &contconf;
         ods output Summary=extreme50u;
      run;

      proc transpose data=extreme5u(drop=NObs &trtvar.) out=extreme5u_t;
      run;

      proc transpose data=extreme95u(drop=NObs &trtvar.) out=extreme95u_t;
      run;

      proc transpose data=extreme50u(drop=NObs &trtvar.) out=extreme50u_t;
      run;

      data percentiles_u;
         attrib perctls_u0 perctls_u1 length=$30.;
         merge extreme5u_t(rename=(col1=p05_0 col2=p05_1))
               extreme95u_t(rename=(col1=p95_0 col2=p95_1))
               extreme50u_t(rename=(col1=p50_0 col2=p50_1));

         if p50_0 > 10 then perctls_u0 = compress(put(p50_0,7.)) || " [" || compress(put(p05_0,7.)) || " - " || compress(put(p95_0,7.)) || "]";
         else perctls_u0 = compress(put(p50_0,7.1)) || " [" || compress(put(p05_0,7.1)) || " - " || compress(put(p95_0,7.1)) || "]";

         if p50_1 > 10 then perctls_u1 = compress(put(p50_1,7.)) || " [" || compress(put(p05_1,7.)) || " - " || compress(put(p95_1,7.)) || "]";
         else perctls_u1 = compress(put(p50_1,7.1)) || " [" || compress(put(p05_1,7.1)) || " - " || compress(put(p95_1,7.1)) || "]";
         keep perctls_u0 perctls_u1;
      run;

      data for_table_unw;
         merge for_table_unw
               percentiles_u;
      run;

      proc append base=stacked_tabunw data=for_table_unw;
      run;

      * Variance Inflation;
      proc means data=sd_data n sum uss  nway noprint;
         class &trtvar;
         var wvar;
         output out=temp_vi n= sum= uss= / autoname;
      run;

      options mergenoby=nowarn;
      data vi;
         merge temp_vi(where=(&trtvar.=0) keep=&trtvar wvar_N wvar_Sum wvar_USS rename=(wvar_N=N0 wvar_Sum=S0 wvar_USS=SS0))
               temp_vi(where=(&trtvar.=1) keep=&trtvar wvar_N wvar_Sum wvar_USS rename=(wvar_N=N1 wvar_Sum=S1 wvar_USS=SS1));

         vi = (1/((1/N1) + (1/N0)))*(SS1/(S1**2) + SS0/(S0**2));
         ss = N1 + N0;

         ess = (1/SS1)*(S1**2) + (1/SS0)*(S0**2);
         ess = round(ess,1);

         subgn = &zz;
         keep subgn vi ss ess;
      run;

      proc means data=sd_data n sum uss  nway noprint;
         var wvar;
         output out=temp_vio sum= uss= / autoname;
      run;

      data temp_vio;
         set temp_vio;
         ess = (1/wvar_USS)*(wvar_Sum**2);
      run;

      proc append base=stacked_vi data=vi;
      run;

      proc append base=stacked_vio data=temp_vio;
      run;
      %end;

   * Add VI/SS to Tables / Index Confounders;
   data stacked_sd;
      merge stacked_sd
            stacked_vi;
      by subgn;

      retain confn;
      if first.subgn then confn=0;
      confn + 1;

   run;

   data stacked_sdunw;
      merge stacked_sdunw
            stacked_vi;
      by subgn;

      retain confn;
      if first.subgn then confn=0;
      confn + 1;

      drop vi ess;
   run;

   * Table;
   data table;
      attrib stats1_u stats0_u stats1_w stats0_w length=$40.;
      merge stacked_tabunw(keep=subgn VarName col1 col2 standiff standifft perctls_u0 perctls_u1 rename=(col1=umean1 col2=umean0 standiff=usd standifft=usdt))
            stacked_tab(keep=col1 col2 standiff standifft perctls_w0 perctls_w1 rename=(col1=wmean1 col2=wmean0 standiff=wsd standifft=wsdt));

      if perctls_u1 ^= "" then do;
         if input(umean1,best8.) > 100 then stats1_u = compress(put(input(umean1,best8.),7.)) || ", " || perctls_u1;
         else stats1_u = compress(put(input(umean1,best8.),7.1)) || ", " || perctls_u1;
         end;
      else if umean1 ^= "" then do;
         if input(umean1,best8.) > 5 then stats1_u = compress(put(input(umean1,best8.),6.1)) || "%";
         else stats1_u = compress(put(input(umean1,best8.),6.2)) || "%";
         end;

      if perctls_u0 ^= "" then do;
         if input(umean0,best8.) > 100 then stats0_u = compress(put(input(umean0,best8.),7.)) || ", " || perctls_u0;
         else stats0_u = compress(put(input(umean0,best8.),7.1)) || ", " || perctls_u0;
         end;
      else if umean0 ^= "" then do;
         if input(umean0,best8.) > 5 then stats0_u = compress(put(input(umean0,best8.),6.1)) || "%";
         else stats0_u = compress(put(input(umean0,best8.),6.2)) || "%";
         end;

      if perctls_w1 ^= "" then do;
         if input(wmean1,best8.) > 100 then stats1_w = compress(put(input(wmean1,best8.),7.)) || ", " || perctls_w1;
         else stats1_w = compress(put(input(wmean1,best8.),7.1)) || ", " || perctls_w1;
         end;
      else if wmean1 ^= "" then do;
         if input(wmean1,best8.) > 5 then stats1_w = compress(put(input(wmean1,best8.),6.1)) || "%";
         else stats1_w = compress(put(input(wmean1,best8.),6.2)) || "%";
         end;

      if perctls_w0 ^= "" then do;
         if input(wmean0,best8.) > 100 then stats0_w = compress(put(input(wmean0,best8.),7.)) || ", " || perctls_w0;
         else stats0_w = compress(put(input(wmean0,best8.),7.1)) || ", " || perctls_w0;
         end;
      else if wmean0 ^= "" then do;
         if input(wmean0,best8.) > 5 then stats0_w = compress(put(input(wmean0,best8.),6.1)) || "%";
         else stats0_w = compress(put(input(wmean0,best8.),6.2)) || "%";
         end;
   run;

   %if %length(&trtfmt) > 0 %then %do;
      proc format library=work cntlout=forms;
         select &trtfmt;
      run;

      data _nothing_;
         set forms;
         if start=1 then call symput("col1",%quote(label));
         else if start=0 then call symput("col0",%quote(label));
      run;
      %end;

   %else %do;
      %let col1 = %str(TRT=1);
      %let col0 = %str(TRT=0);
      %end;

   options orientation=landscape;
   ods rtf file="&outdir.&name._table.rtf";
   ods escapechar='!';
   proc report data=table nowd spanrows split = '|' missing
      style(report)={just=center}
      style(lines)=header{background=white font_size=9pt font_face="Arial" font_weight=bold just=left}
      style(header)=header{background=white font_size=9pt font_face="Arial" font_weight=medium}
      style(column)=header{background=white font_size=9pt font_face="Arial" font_weight=medium protectspecialchars = off};

      columns subgn VarName ("Unweighted" stats1_u stats0_u usdt)
                            ("Weighted"   stats1_w stats0_w wsdt);

      define subgn / display "Subgroup" style(column)={just=left cellwidth=10%} order order=data;
      define VarName / display "Confounder" style(column)={just=left cellwidth=17%};

      define stats1_u /display "&col1 | Mean, Median [P!{sub 5} - P!{sub 95}]" style(column)=[just=center cellwidth=15%];
      define stats0_u /display "&col0 | Mean, Median [P!{sub 5} - P!{sub 95}]" style(column)=[just=center cellwidth=15%];
      define usdt /display "Stand.| Diff." style(column)=[just=center cellwidth=6%];
      define stats1_w /display "&col1 | Mean, Median [P!{sub 5} - P!{sub 95}]" style(column)=[just=center cellwidth=15%];
      define stats0_w /display "&col0 | Mean, Median [P!{sub 5} - P!{sub 95}]" style(column)=[just=center cellwidth=15%];
      define wsdt /display "Stand.| Diff." style(column)=[just=center cellwidth=6%];
      %if %length(&conffmt) > 0 %then %do;
         format VarName &conffmt..;
         %end;

      %if %length(&subgfmt) > 0 %then %do;
         format subgn &subgfmt..;
         %end;
      %else %do;
         format subgn subgbf.;
         %end;
   run;
   ods rtf close;

   * ConnectS plot;
   %connectS(dsp=stacked_sdunw,stdvble=standiff_unw,confidx=VarName,subgidx=subgn,
             gpath=&outdir.,imgname=&name._connect_unweighted
             %if %length(&conffmt) > 0 %then %do;
                ,conff=&conffmt
                %end;
             %if %length(&subgfmt) > 0 %then %do;
                ,subgf=&subgfmt
                %end;
             %else %do;
                ,subgf=subgbf
                %end;
             %if &ss = Y %then %do;
                ,sampsize=ss
                %end;
             %if %length(&hgt) > 0 %then %do;
                ,h=&hgt
                %end;
             %if %length(&wdt) > 0 %then %do;
                ,w=&wdt
                %end;
             %if %length(&rf) > 0 %then %do;
                ,refs=&reflines
                %end;
             %if %length(&fsize) > 0 %then %do;
                ,fs=&fsize
                %end;
             %if %length(&csize) > 0 %then %do;
                ,cs=&csize
                %end;
             %if %length(&color) > 0 %then %do;
                ,col=&color
                %end;
             );

   %connectS(dsp=stacked_sd,stdvble=standiff,confidx=VarName,subgidx=subgn,
             gpath=&outdir.,imgname=&name._connect_weighted
             %if %length(&conffmt) > 0 %then %do;
                ,conff=&conffmt
                %end;
             %if %length(&subgfmt) > 0 %then %do;
                ,subgf=&subgfmt
                %end;
             %else %do;
                ,subgf=subgbf
                %end;
             %if &ss = Y %then %do;
                ,sampsize=ess
                %end;
             %if &vi = Y %then %do;
                ,vivble=vi
                %end;
             %if %length(&hgt) > 0 %then %do;
                ,h=&hgt
                %end;
             %if %length(&wdt) > 0 %then %do;
                ,w=&wdt
                %end;
             %if %length(&rf) > 0 %then %do;
                ,refs=&reflines
                %end;
             %if %length(&fsize) > 0 %then %do;
                ,fs=&fsize
                %end;
             %if %length(&csize) > 0 %then %do;
                ,cs=&csize
                %end;
             %if %length(&color) > 0 %then %do;
                ,col=&color
                %end;
             );

   * Histograms;
   %if %length(&propens) > 0 %then %do;
      %histolap(dsq=stacked,trt=&trtvar,subgidx=subgn,xvar=&propens,wvble=wvar,
                gpath=&outdir.,imgname=&name._propensity_histograms,xmi=0,xma=1,xst=0.2,
                xlabel=%str(Weighted Propensity)
                %if %length(&subgfmt) > 0 %then %do;
                   ,subgf=&subgfmt
                   %end;
                %else %do;
                   ,subgf=subgbf
                   %end;
                %if %length(&trtfmt) > 0 %then %do;
                   ,trtf=&trtfmt
                   %end;
                %if %length(&color) > 0 %then %do;
                   ,col=&color
                   %end;
                );
      %end;

   %if %length(&prognos.) > 0 %then %do;
      proc sql noprint;
         select max(&prognos.) into :maxprog from stacked;
         select min(&prognos.) into :minprog from stacked;
      quit;

      %if %length(&hxlabel) > 0 %then %do;
         %let xla = &hxlabel;
         %end;

      %else %do;
         %let xla = %str(&prognos.);
         %end;

      %histolap(dsq=stacked,trt=&trtvar,subgidx=subgn,xvar=&prognos,wvble=wvar,
                gpath=&outdir.,imgname=&name._&prognos._histograms,xlabel=&xla
                %if %length(&subgfmt) > 0 %then %do;
                   ,subgf=&subgfmt
                   %end;
                %else %do;
                   ,subgf=subgbf
                   %end;
                %if %length(&trtfmt) > 0 %then %do;
                   ,trtf=&trtfmt
                   %end;
                %if %length(&color) > 0 %then %do;
                   ,col=&color
                   %end;
                %if %length(&xmin) > 0 %then %do;
                   ,xmi=&xmin
                   %end;
                %else %do;
                   ,xmi=&minprog
                   %end;
                %if %length(&xmax) > 0 %then %do;
                   ,xma=&xmax
                   %end;
                %else %do;
                   ,xma=&maxprog
                   %end;
                %if %length(&xstep) > 0 %then %do;
                   ,xst=&xstep
                   %end;
                );
      %end;
   ods listing;
   proc datasets library=work;
      delete _all_;
   quit;
%mend diagSGA;

* ######################################################################################### *
* ################################### stddiff macro ####################################### *
* ######################################################################################### *;
* modified STDDIFF macro;
%macro stddiff(inds=, groupvar =, numvars=, charvars=, wtvar=, stdfmt=8.4, outds = stddiff_result);
   * create a table to store stddiff and means/percentages;
   proc sql;
      create table &outds.
         (VarName char(32),
   		Stddiff char (10),
			col1 char(10),
			col2 char(10)
         );
   quit;

   * delete records if the group variable is missing;
   data base_data;
 	   set &inds.;
 	   where &GroupVar. ne .;
   run;

   * remove leading or tailing blanks;
   %let groupvar = %sysfunc(strip(&GroupVar.));

   * part 1: compare continuous variables *;
   %if %length(&numvars.) > 0 %then %do;

      * remove multiple blanks and get the total number of continuous variables;
	   %let numvar = %sysfunc(compbl(&numvars.));
	   %let numvar = %sysfunc(strip(&numvar.));
	   %let n_convar = %sysfunc(countc(&numvar.,' '));
	   %let n_convar = %eval(&n_convar. + 1);

      * summarize variables one-by-one;
	   %do ii = 1 %to &n_convar.;

    	   %let convar = %sysfunc(scan(&numvar.,&ii.,' '));

		   %if %index(&convar., /r) > 0 %then %do;
    		   %let convar = %sysfunc(scan(&convar.,1,'/'));
    		   %let convar = %sysfunc(strip(&convar.));

    		   data temp_1;
     			   set base_data (keep = &groupvar. &convar. &wtvar.);
    		   run;

            * rank a variable;
    		   proc rank data=temp_1 out=temp_2;
         	   var &convar.;
        		   ranks rank_&convar.;
    		   run;

            * get ranked-mean and sd;
			   proc means data = temp_2 vardef=weight;
				   class &groupvar.;
				   var rank_&convar.;
				   weight &wtvar.;
				   output out = temp_3 mean = _mean_  std = _std_;
			   run;

			   data temp_3;
				   set temp_3;
				   where _type_ = 1;
			   run;

			   proc sort data = temp_3;
				   by &groupvar.;
			   run;
	     	   %end;

		   %else %do;
	    	   %let convar = %sysfunc(strip(&convar.));
	    	   data temp_1;
	    	 	   set base_data (keep = &groupvar. &convar. &wtvar.);
	    	   run;

	    	   data temp_2;
	     		   set temp_1;
	    	   run;

	         * get mean and sd;
			   proc means data = temp_2 vardef=weight;
				   class &groupvar.;
				   var &convar.;
				   weight &wtvar.;
				   output out = temp_3 mean = _mean_  std = _std_;
			   run;

			   data temp_3;
				   set temp_3;
				   where _type_ = 1;
			   run;

			   proc sort data = temp_3;
				   by &groupvar.;
			   run;
    		   %end;

         * calculate stddiff *;
	      proc sql;
	    	   create table temp_4 as
	    		   select (a._mean_ - b._mean_)/sqrt((a._std_**2 + b._std_**2)/2) as d,
				           a._mean_ as e,
				           b._mean_ as f
	    		   from temp_3(where = (&groupvar = 1)) as a,
	      		     temp_3(where = (&groupvar = 0)) as b;
	      quit;

	      data temp_5;
	   	   set temp_4;
	         stddiff = compress(put(d,&stdfmt.));
			   col1v = compress(put(e,5.1));
			   col2v = compress(put(f,5.1));
	         keep col1v col2v stddiff;
	      run;

	      * insert into std table;
	      proc sql noprint;
	    	   select stddiff into: std_value from temp_5;
			   select col1v into: col1val from temp_5;
			   select col2v into: col2val from temp_5;
	    	   insert into &outds.  values("&convar.","&std_value.","&col1val","&col2val");
	      quit;

	      * delete temporary data sets *;
	      proc datasets lib = work nodetails nolist;
	         delete temp_1 - temp_5;
	      quit;
   	   %end;
   %end;

* part 2: compare categorical variables     *;
%if %length(&charvars.) > 0 %then %do;
	%let n_charvar = %sysfunc(countw(&charvars.));

   * get column percents for each levels of the variable by the group;
	%do jj = 1 %to &n_charvar.;
   		%let char_var = %scan(&charvars., &jj.);
   		%let char_var = %sysfunc(strip(&char_var.));
		  data temp_1;
		   	set base_data (keep = &groupvar. &char_var. &wtvar.);
		  run;

		  proc sql;
		   	create table temp_2 as
		   	select distinct &char_var. as &char_var.
		   	from temp_1
			where &char_var. is not missing;
		  quit;

		  proc sql noprint;
		   	select count(*) into :_mylevel_ from temp_2;
		  quit;

   	  %let _mylevel_ = %sysfunc(strip(&_mylevel_.));

		  data temp_3;
		   	set temp_2;
		   		do &groupvar. = 0,1 ;
		    	      output;
		   		   end;
		  run;

		  ods output CrossTabFreqs = temp_4;
  		  proc freq data = temp_1;
   	     table &char_var. * &groupvar.;
			     %if %length(&wtvar.) > 0 %then %do;
				  weight &wtvar.;
				  %end;
  		  run;

	  	  proc sql;
	        create table  temp_5 as select a.*, b.ColPercent
	   		from temp_3 as a
	   		left join temp_4 as b
	   		on a.&groupvar. = b.&groupvar. and a.&char_var. = b.&char_var.;
	  	  quit;

	  	  data temp_6;
	        set temp_5;
	   	  if ColPercent = . then ColPercent = 0;
	  	  run;

  		  proc sort data = temp_6 out = catfreq;
   		  by &groupvar. &char_var.;
  		  run;

  		  proc datasets lib = work nodetails nolist;
   		  delete  temp_1 - temp_6;
  		  quit;

        * if a categorical variable only has one level: 0 or 1;
        * stddiff = 0 *;
		  %if &_mylevel_. = 1 %then %do;
  		     proc sql noprint;
   		     insert into &outds. values("&char_var.", "0"," "," ");
  			  quit;
   		  %end;

		  %else %if &_mylevel_. = 2 %then %do;

  		  data temp_7;
   	     set catfreq;
  			  where &char_var. = 1;
   		  ColPercent = ColPercent/100;
  		  run;

  		  proc sql;
           create table temp_8 as
   		  select (a.ColPercent - b.ColPercent)/(sqrt((a.ColPercent*(1-a.ColPercent) +
     					 b.ColPercent*(1-b.ColPercent))/2)) as d,
				       a.ColPercent as e,
				       b.ColPercent as f
   				    from temp_7(where = (&groupvar = 1)) as a,
     		    	    temp_7(where = (&groupvar = 0)) as b;
 		  quit;

  		  data temp_9;
           set temp_8;
           stddiff = compress(put(d,&stdfmt.));
  			  col1v = compress(put(100*e,5.1));
			  col2v = compress(put(100*f,5.1));
	        keep col1v col2v stddiff;
        run;

  		  proc sql noprint;
           select stddiff into: std_value from temp_9;
			  select col1v into: col1val from temp_9;
			  select col2v into: col2val from temp_9;
	    	  insert into &outds.  values("&char_var.","&std_value.","&col1val","&col2val");
        quit;

  		  proc datasets lib = work nodetails nolist;
   	     delete  temp_7 temp_8 temp_9;
  		  quit;
    	  %end;

		%else %if &_mylevel_. > 2 %then %do;
   	   %let _k_ = %eval(&_mylevel_. - 1);
   		%let _k_ = %sysfunc(strip(&_k_.));

  			data temp_7;
   		   set catfreq;
  				by &groupvar.;
   			if last.&groupvar. then delete;
   			ColPercent = ColPercent/100;
  			run;

  			proc sql noprint;
   		   select ColPercent format=9.6 into :tlist separated by ' '
				from temp_7 where &groupvar. = 1;
  				select ColPercent format=9.6 into :clist separated by ' '
				from temp_7 where &groupvar. = 0;
  			quit;

      * vector T, C and T-C *;
  			data t_1;
   		   array t{*}  t1- t&_k_.   (&tlist.);
   			array c{*}  c1- c&_k_.   (&clist.);
   			array tc{*} tc1 - tc&_k_. ;
   			do i = 1 to dim(t);
    			   tc{i} = t{i} - c{i};
   				end;
   			drop i;
  			run;

      * each column has one element of a S covariance matrix (k x k) *;
			%let _dm = ;
			%let _dm = %eval(&_k_.*&_k_.);
  			data covdata;
   		   array t{*}  t1- t&_k_.  (&tlist.);
   			array c{*}  c1- c&_k_.   (&clist.);
   			array cv{&_k_.,&_k_.} x1 -x&_dm.;
   			do i = 1 to &_k_.;
    				do j = 1 to &_k_.;
     					if i = j then do;
      				   cv{i,j} = 0.5*(t{i}*(1-t{i}) + c{i}*(1-c{i}));
      					end;
     					else do;
      				   cv{i,j} = -0.5 * (t[i] * t[j] + c[i] * c[j]);
      					end;
    					if cv{&_k_.,&_k_.] ne . then output;
    				   end;
  				   end;
  			 run;

  			 proc transpose data = covdata(keep = x1 -x&_dm.) out = covdata_1;
  			 run;

  			 data covdata_2;
   		    set covdata_1;
   			 retain id gp 1;
   			 if mod(_n_ - 1,&_k_.) = 0 then gp = gp + 1;
  			 run;

		  	 proc sort data = covdata_2 ;
		       by gp id;
		  	 run;

			 data covdata_3;
		       set covdata_2;
		   	 by gp id;
		   	 retain lp;
		   	 if first.gp then lp = 0;
		   	 lp = lp+1;
		  	 run;

          * transpose to a S variance-covariance matrix format *;

		  	 data covdata_4;
		       set covdata_3;
		   	 retain y1-y&_k_.;
		   	 array cy{1:&_k_.} y1-y&_k_.;
		   	 by gp id;
		   	 if first.gp then do;
		    	    do k = 1 to &_k_.;
		     		    cy{k} = .;
				       end;
		   		 end;
		   	 cy{lp} = col1;
		   	 if last.gp then output;
		   	 keep y:;
		  	 run;

          * get inverse of S matrix *;
		    data A_1;
		       set covdata_4;
		       array _I{*} I1-I&_k_.;
		       do j=1 to &_k_.;
		          if j=_n_ then _I[j]=1;
		          else _I[j]=0;
		          end;
		       drop j;
		    run;

          * solve the inverse of the matrix *;
          %macro inv;
    	       %do j=1 %to &_k_.;
    		       proc orthoreg data=A_1 outest=A_inv_&j.(keep=y1-y&_k_.) noprint singular=1E-16;
     			       model I&j=y1-y&_k_. /noint;
    		       run;
    		       quit;
    	          %end;

   		    data A_inverse;
    		       set %do j=1 %to &_k_.;
     		          A_inv_&j
     	             %end;;
   		    run;
          %mend;
          %inv;

  		    proc transpose data=A_inverse out=A_inverse_t;
  		    run;

          * calculate the mahalanobis distance *;
  		    data t_2;
   		    set A_inverse_t;
   			 array t{*}  t1- t&_k_.  (&tlist.);
   			 array c{*}  c1- c&_k_.  (&clist.);
   			 i = _n_;
   			 trt = t{i};
   			 ctl = c{i};
   			 tc = t{i} - c{i};
  		    run;

		    data t_3;
   		    set t_2;
   			 array aa{&_k_.} col1 - col&_k_.;
   			 array bb{&_k_.} bb1- bb&_k_.;
   			 do i = 1 to &_k_.;
    			    bb{i} = aa{i}*tc;
   			    end;
  		    run;

  		    proc summary data = t_3 ;
   		    var bb1-bb&_k_.;
   			 output out = t_4 sum =;
  		    run;

  		    data t_5;
   		    merge t_1 t_4;
   			 array d1{*} tc1- tc&_k_. ;
   			 array d2{*} bb1-bb&_k_.;
   			 array d3{*} y1-y&_k_.;
   			 do i = 1 to &_k_.;
   			    d3{i} = d1{i}*d2{i};
   			    end;
   			 d = sqrt(sum(of y1-y&_k_.));
   			 stddiff = compress(put(d,&stdfmt.));

   			 keep stddiff;
  		    run;

  		    proc sql noprint;
   		    select  stddiff into: std_value from t_5;
   			 insert into &outds.  values("&char_var.","&std_value."," "," ");
  		    quit;

		    data t_6;
			    merge catfreq(where=(&groupvar=1) rename=ColPercent=col1v)
			          catfreq(where=(&groupvar=0) rename=ColPercent=col2v);

			    keep &char_var col1v col2v;
          run;

		    data t_7;
			    attrib VarName length=$32.;
			    set t_6;

 			    col1 = compress(put(col1v,5.1));
			    col2 = compress(put(col2v,5.1));
			    VarName = &char_var;

             keep VarName col1 col2;
          run;

		    data &outds.;
			    set &outds
			        t_7;
          run;

  		    proc datasets lib = work nodetails nolist;
   		    delete covdata covdata_1 covdata_2 covdata_3 covdata_4
      		        A_1 A_inverse A_inverse_t t_1 t_2 t_3 t_4 t_5 t_6 t_7 A_inv_:;
  		    quit;
	       %end;
	   %end;
   %end;

   proc datasets lib = work nodetails nolist;
     delete Catfreq  Base_data temp_7;
   quit;
   title;
%mend;


* ######################################################################################### *
* ################################## connectS macro ####################################### *
* ######################################################################################### *;
data myrattrmap;
   retain id "myid";
   length min $ 5 max $ 5;
   input min $ max $ colormodel1 $ colormodel2 $ altcolormodel1 $ altcolormodel2 $;
   datalines;
   0 0.10 green yellow green yellow
0.10 0.25 yellow red yellow red
0.25 0.50 red darkred red darkred
;
run;

%macro connectS(dsp=,stdvble=,confidx=,subgidx=,vivble=,sampsize=,gpath=,
                imgname=connectS_plot,h=,w=,tit=,refs=,subgf=,conff=,fs=8,cs=,col=);

%if &dsp=stacked_sdunw %then %do;
   %let sslab = "Sample Size";
   %end;

%else %do;
   %let sslab = "Effective Sample Size";
   %end;

data ds1;
   set &dsp;

   if abs(&stdvble) > 0.20 then stdcat=3;
   else if 0.15 < abs(&stdvble) <= 0.20 then stdcat=2;
   else if 0.1 < abs(&stdvble) <= 0.15 then stdcat=1;
   else if abs(&stdvble) <= 0.1 then stdcat=0;

   sd_abs = abs(&stdvble);
run;

proc sort data=ds1;
   by &confidx &subgidx;
run;

data order;
   do stdcat = 0 to 3;
      output;
   end;
run;

data to_plot;
   set order
       ds1;
run;

proc sql noprint;
   select max(&subgidx) into :nsubgs from to_plot;
quit;

proc sort data=to_plot out=countconf nodupkey;
   by &confidx;
run;

proc sql noprint;
   select n(&confidx) into :nconf from countconf;
quit;

proc format;
   value catlev 0="ASMD (*ESC*){unicode '2264'x} 0.10" 1="ASMD 0.10 - 0.15" 2="ASMD 0.15 - 0.20" 3="ASMD > 0.20";
run;

%let fss = %sysevalf(0.75*&fs);

%if %length(&w) = 0 %then %do;
   %let ww = %sysevalf(0.4*&nconf);
   %end;

%else %if %length(&w) > 0 %then %do;
   %let ww = &w;
   %end;

%if %length(&h) = 0 %then %do;
   %let hh = %sysevalf(0.4*&nsubgs);
   %end;

%else %if %length(&h) > 0 %then %do;
   %let hh = &h;
   %end;

%if %length(&cs) = 0 %then %do;
   %let ccs = 12;
   %end;

%else %if %length(&cs) > 0 %then %do;
   %let ccs = &cs;
   %end;

   ods listing image_dpi=600 gpath="&gpath";
   ods graphics on / imagename="&imgname" noborder height=&hh.in width=&ww.in;

   proc sgplot data=to_plot;
      %if %length(&tit) > 0 %then %do;
         title "&tit";
         %end;

      styleattrs datacolors=(white grayC0 gray70 black);
      scatter y=&subgidx x=&confidx / markerattrs=(size=&ccs symbol=circlefilled) filledoutlinedmarkers
                                      markeroutlineattrs=(color=black thickness=0.1)
							                 group=stdcat name="plot" grouporder=ascending;

      %if %length(&refs) > 0 %then %do;
         refline &refs / axis=y;
         %end;

      xaxis label="Confounder" valueattrs=(size=&fs.pt) fitpolicy=rotate type=discrete;
      yaxis reverse values=(1 to &nsubgs by 1) label=%str("Subgroups") valueattrs=(size=&fs.pt);

      %if %length(&sampsize) > 0 %then %do;
         yaxistable &sampsize. / y=&subgidx label=&sslab labelattrs=(size=&fss.pt) location=inside position=right valueattrs=(size=&fs.pt) valuejustify=center pad=(left=0px) stat=mean;
         %end;

      %if %length(&vivble) > 0 %then %do;
         yaxistable &vivble / y=&subgidx label="Variance Inflation" labelattrs=(size=&fss.pt) location=inside position=right valueattrs=(size=&fs.pt) valuejustify=center pad=(left=0px) stat=mean;
         %end;

      keylegend  "plot" / noborder valueattrs=(size=&fs.pt) title="Absolute Standardized Mean Difference" titleattrs=(size=&fs.pt) down=1;
      format stdcat catlev.;

      %if %length(&subgf) > 0 %then %do;
         format &subgidx &subgf..;
         %end;

      %if %length(&conff) > 0 %then %do;
         format &confidx &conff..;
         %end;

      %if %length(&vivble) > 0 %then %do;
         format &vivble 6.2;
         %end;
   run;

%if &col = Y %then %do;
   ods listing image_dpi=600 gpath="&gpath";
   ods graphics on / imagename="&imgname._color" noborder height=&hh.in width=&ww.in;

   proc sgplot data=to_plot rattrmap=myrattrmap;
      %if %length(&tit) > 0 %then %do;
         title "&tit";
         %end;

      scatter y=&subgidx x=&confidx / markerattrs=(size=&ccs symbol=circlefilled) filledoutlinedmarkers
                                      markeroutlineattrs=(color=white thickness=0.1)
							                 colorresponse=sd_abs rattrid=myid name="plot";

      %if %length(&refs) > 0 %then %do;
         refline &refs / axis=y;
         %end;

      xaxis label="Confounder" valueattrs=(size=&fs.pt) fitpolicy=rotate type=discrete;
      yaxis reverse values=(1 to &nsubgs by 1) label=%str("Subgroups") valueattrs=(size=&fs.pt);

      %if %length(&sampsize) > 0 %then %do;
         yaxistable &sampsize. / y=&subgidx label=&sslab labelattrs=(size=&fss.pt) location=inside position=right valueattrs=(size=&fs.pt) valuejustify=center pad=(left=0px) stat=mean;
         %end;

      %if %length(&vivble) > 0 %then %do;
         yaxistable &vivble / y=&subgidx label="Variance Inflation" labelattrs=(size=&fss.pt) location=inside position=right valueattrs=(size=&fs.pt) valuejustify=center pad=(left=0px) stat=mean;
         %end;

      gradlegend "plot" / noborder title="Absolute Standardized Mean Difference" titleattrs=(size=&fs.pt) position=bottom;

      %if %length(&subgf) > 0 %then %do;
         format &subgidx &subgf..;
         %end;

      %if %length(&conff) > 0 %then %do;
         format &confidx &conff..;
         %end;

      %if %length(&vivble) > 0 %then %do;
         format &vivble 5.1;
         %end;
   run;
   %end;
%mend connectS;

* ######################################################################################### *
* ################################## histolap macro ####################################### *
* ######################################################################################### *;
%macro histolap(dsq=,trt=,subgidx=,xvar=,wvble=,gpath=,imgname=,col=N,rows=,cols=,tit=,subgf=,trtf=,
                xmi=,xma=,xst=,xlabel=%str(Weighted Propensity));
data histo;
   set &dsq;
run;

proc sort data=histo;
   by &subgidx &trt;
run;

proc means data=histo noprint nway;
   var &subgidx;
   output out=count max(&subgidx)=nsubgs;
run;

data count;
   set count;

   columns_n = ceil(sqrt(nsubgs));
   rows_n = ceil(nsubgs/columns_n);
   call symput('rowsx', compress(put(rows_n, best9.)));
   call symput('colsx', compress(put(columns_n, best9.)));
run;

%if %length(&rows) = 0 or %length(&cols) = 0 %then %do;
   %let rows = &rowsx;
   %let cols = &colsx;
   %end;

%let ww = %sysevalf(1.5*&cols);
%let hh = %sysevalf(1.5*&rows);

%if &col = N %then %do;
   %let co_scale = (gray40 gray90);
   %end;

%else %do;
   %let co_scale = ("#D53333" "#6B7AB9");
   %end;

ods graphics on / imagename="&imgname" noborder height=&hh.in width=&ww.in;
ods listing image_dpi=300 gpath="&gpath";

proc sgpanel data=histo;

   %if %length(&tit) > 0 %then %do;
      title "&tit";
      %end;

   styleattrs datacolors=&co_scale;
   panelby &subgidx / columns=&cols rows=&rows novarname spacing=5 headerattrs=(size=9) noborder
                      noheaderborder skipemptycells uniscale=column;
   histogram &xvar / group=&trt transparency=0.6 name="hist" weight=&wvble nooutline;
   %if %length(&xst) > 0 %then %do;
      colaxis values=(&xmi to &xma by &xst) label="&xlabel" labelattrs=(weight=normal) valueattrs=(size=7);
      %end;
   %else %do;
      colaxis values=(&xmi to &xma) label="&xlabel" labelattrs=(weight=normal) valueattrs=(size=7);
      %end;
   rowaxis labelattrs=(weight=normal);

   %if %length(&subgf) > 0 %then %do;
      format &subgidx &subgf..;
      %end;

   %if %length(&trtf) > 0 %then %do;
      format &trt &trtf..;
      %end;

   keylegend "hist" / noborder valueattrs=(size=10) outerpad=(top=0.15in) autooutline fillheight=10;
run;
%mend histolap;
