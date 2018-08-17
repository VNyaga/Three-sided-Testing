/*
CREATED:	25 Jul 2018
AUTHOR:		Victoria N Nyaga
PURPOSE: 	Three-sided hypothesis testing.
VERSION: 	1.0.0
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
UPDATES
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
DATE:						DETAILS:

*/


cap program drop threesidedtest
program define threesidedtest, rclass
version 14.0
	#delimit ;
	syntax anything(name=data id="data"), delta_l(string) delta_u(string) 
	[noKEY KEYOptions(string asis) noSTAT vline(string) ADDploti(string asis) noGRaph BColors(string) 
	alpha(real 0.05) dp(integer 3) 
	*
	]
	;
	#delimit cr
	
	foreach num of local data {
		cap confirm integer number `num'
		if _rc != 0 {
			di as error "`num' found where integer expected"
			exit
		}
	}
	
	local len: word count `data'
	if `len' != 4 {
		di as error "Specify full data: a b c d"
		exit
	}
	
	tokenize `data'
	local a = `1'
	local b = `2'
	local c = `3'
	local d = `4'
	
	local n = `a' + `b' + `c' + `d'
	
	local p1 = (`a' + `b')/`n'
	local p0 = (`a' + `c')/`n'
	local p01 = `c'/`n'
	local p10 = `b'/`n'
	local p00 = `d'/`n'
	local RR = `p1'/`p0'
	local level = (1 - `alpha')*100
	
	
	/*Test statistics*/
	zvalue `p1' `p0' `p10' `p01'  `p00' `n', delta(`delta_l') 
	local z_l = r(z)	
	
	zvalue `p1' `p0' `p10' `p01'  `p00' `n', delta(`delta_u') 
	local z_u = r(z)
	
	/*P-values*/

	local noninf_p :di %10.`dp'f 1-normprob(`z_l')
	local super_p :di %10.`dp'f 1-normprob(`z_u')

	local equiv_p1 = 1-normprob(`z_l')
	local equiv_p2 = normprob(`z_u')
	local equiv_p = max(`equiv_p1', `equiv_p2')
	
	/*Wald CI
	local inits_l = `RR'*exp(-1*invnorm(1 - `alpha'/2)*sqrt((`b' + `c')/((`a' + `b')*(`a' + `c'))))
	local inits_u = `RR'*exp(invnorm(1 - `alpha'/2)*sqrt((`b' + `c')/((`a' + `b')*(`a' + `c'))))
	*/
	
	mata: cml_ci((`a', `b', `c', `d'), `alpha')
	
	/*Three-sided CI*/
	/*Lower*/
    if (ci[1,1] < `delta_l') {
	
        local newl :di %10.`dp'f ci[1,1]
		local newl =  "(" + strtrim("`newl'")  
		
    } 
	else if ((ci[2,1] <= `delta_l') & (`delta_l' <= ci[1, 1])){
	
        local newl :di %10.`dp'f `delta_l'
		local newl =  "[" + strtrim("`newl'") 
        
    } 
	else if((`delta_l' < ci[2,1]) & (ci[2,1] < `delta_u')){
        
        local newl :di %10.`dp'f ci[2,1]
		local newl =  "(" + strtrim("`newl'")  
        
    } 
	else if (`delta_u' <= ci[2,1]){
        
        local newl :di %10.`dp'f `delta_u'
		local newl =  "(" + strtrim("`newl'") 
        
    }
    
	/*Upper*/
    if (ci[2,2] <= `delta_l'){
		
		local newu :di %10.`dp'f `delta_l'
		local newu =  strtrim("`newu'") + ")" 
        
    } 
	else if((`delta_l' < ci[2, 2]) & (ci[2, 2] < `delta_u')){
        
		local newu :di %10.`dp'f ci[2,2]
		local newu =  strtrim("`newu'") + ")"
        
    } 
	else if((ci[1, 2] <= `delta_u') & (`delta_u' <= ci[2,2])){
        
		local newu :di %10.`dp'f `delta_u'
		local newu =  strtrim("`newu'") + "]"
        
    } 
	else if(`delta_u' < ci[1, 2]){
        
		local newu :di %10.`dp'f ci[1, 2]
		local newu =  strtrim("`newu'") + ")"
    }
	
	local ci = "`newl'" + ", " + "`newu'"

	/*Display*/
	mat mydata =(`a', `b', `a' + `b')\(`c', `d', `c' + `d')\(`a' + `c', `b' + `d', `n')
	mat colnames mydata = "Controls: Positive" "Controls: Negative" "Controls: Total"
	mat rownames mydata = Positive Negative Total
	/*Data*/
	di
	#delimit ;
	matlist mydata,  rowtitle(Cases) ///
				cspec(& %8s |  %10.0f &  %10.0f |  %10.0f &) ///
				rspec(&-&-&) nohalf
	;
	#delimit cr
	
	/*RR and CI*/
	di
	di as txt "Proportion with factor"
	di as txt _skip(13) "Cases"  _skip(3) as res %10.`dp'f `p1'
	di as txt _skip(13) "Controls" as res %10.`dp'f `p0' _skip(3) as txt "`level'% Conf. Interval"
	di as txt _skip(13) "{hline 17}" _skip(3) "{hline 17}"
	di as txt _skip(13) "Ratio" _skip(3) as res %10.`dp'f `RR' _skip(3) as res "`newl'" _skip(3) "`newu'"
	
	
	/*Equivalence*/
	local z_lf: di %10.`dp'f `z_l'
	local z_lf = strtrim("`z_lf'")
	
	local z_uf: di %10.`dp'f `z_u'
	local z_uf = strtrim("`z_uf'")
	
	di
	di as txt "H0: Not similar (RR < `delta_l' or RR > `delta_u') vs H1: Similar (`delta_l' <= RR <= `delta_u')"
	di as txt _skip(13) "Test statistic(1) = " as res  `z_lf'
	di as txt _skip(13) "Test statistic(2) = " as res  `z_uf'
	if (`equiv_p') < `=0.1^`dp'' {
		local equiv_p =	"`=0.1^`dp''"
		di as txt _skip(13) "P-value" _skip(11) "< " as res `equiv_p'
	}
	else {
		local equiv_p :di %10.`dp'f `equiv_p'
		di as txt _skip(13) "P-value" _skip(11) "= " as res `equiv_p'
	} 
	
	/*Non-Inferiority*/
	di
	di as txt "H0: Inferior (RR < `delta_l') vs H1: Non-inferior(RR >= `delta_l')"
	di as txt _skip(13) "Test statistic " _skip(3) "= " as res `z_lf' 
	if (`noninf_p') < `=0.1^`dp'' {
		local noninf_p =	"`=0.1^`dp''"
		di as txt _skip(13) "P-value" _skip(11) "< " as res `noninf_p'
	}
	else {
		local noninf_p :di %10.`dp'f `noninf_p'
		di as txt _skip(13) "P-value" _skip(11) "= "  as res `noninf_p'
	}

	/*Superiority*/
	di
	di as txt "H0: Non-superior (RR < `delta_u') vs H1: Superior (RR >= `delta_u')"
	di as txt _skip(13) "Test statistic " _skip(3) "= " as res `z_uf'
	if (`super_p') < `=0.1^`dp'' {
		local super_p =	"`=0.1^`dp''"
		di as txt _skip(13) "P-value" _skip(11) "< " as res `super_p'
	}
	else {
		local super_p :di %10.`dp'f `super_p'
		di as txt _skip(13) "P-value" _skip(11) "= "  as res `super_p'
	}
	
	/*Conclusion*/
	di 
	di as txt "Overall conclusion :"
	local syml = substr("`newl'", 1, 1)
	local newl = substr("`newl'", 2, .)
	
	
	local len = strlen("`newu'") - 1
	local symu = substr("`newu'", `=`len' +1', `=`len' +1')
	local newu = substr("`newu'", 1, `len')
	
	if(`newu'  < `delta_l') {
		di as res "New test is inferior" 
	}
	
    if(`newl' > `delta_u') {	
		di as res "New test is superior"
	}
	
    if(( `newl' >= `delta_l') & (`newu' <= `delta_u')) {
		di as res "New test is similar"
	}
	
    if((`newu' < `delta_u') & (`newl' < `delta_l')) {
		di as res "New test is non-superior" 
	}
	
    if((`newu' > `delta_u') & (`newl' > `delta_l')) {
		di as res "New test is non-inferior"
	}
	
	/*Graph*/
	if "`graph'" == "" {
		preserve
		qui {
			clear
			set obs 3
			tempname x
			gen `x' = _n
			
			local delta_uorig = "`delta_u'"
			local delta_u = `=`delta_u''
			local limitl = max(`delta_l' - 0.5, 0)
			local limitu = `delta_u' + 0.5
			
			local midp1 = 0.5*(`limitl' + `delta_l')
			local midp2 = 0.5*(`delta_l' + `delta_u')
			local midp3 = 0.5*(`delta_u' + `limitu')
			
			tempname delta_ly delta_uy limituy
			
			gen `delta_ly' = `delta_l'
			gen `delta_uy' = `delta_u'
			gen `limituy' = `limitu'
		}
		local defaultcolors = "red orange green"
		forvalues i=1/3 {
			local color`i': word `i' of `bcolors'
			if "`color`i''" == "" {
				local color`i': word `i' of `defaultcolors'
			}
		}
				
		local myopts `""yscale(noextend noline)""' 
		local myopts `"`myopts' `"ylabel("", noticks nolabels gextend nogrid)"'"'
		local myopts `"`myopts' `"legend(off)"'"'
		local myopts `"`myopts' `"graphregion(margin(r=0 l=0) color(white))"'"'
		local myopts `"`myopts' `"xlabel(`delta_l' "`delta_l'" 1 "1"  `delta_u' "`delta_uorig'", labc(black) tlcolor(black) axis(1))"'"'
		local myopts `"`myopts' `"xscale(axis(1))"'"'
		local myopts `"`myopts' `"title("Three-sided testing")"'"'
		local myopts `"`myopts' `"xtitle("Ratio")"'"'
		
		foreach opt of local myopts {
			if strpos(`"`options'"', substr(`"`opt'"', 1, 4)) == 0 {
				local options `"`options' `opt'"'
			}
		}
				
		if "`vline'" != "" {
			local vline = "(pci 2.8 `vline' 1 `vline', lcolor(black) lpattern(-))"
		}
		if "`stats'" == "" {
			local RR: di %10.`dp'f `RR'
			local RR = strtrim("`RR'")
			local showstat = `"(scatteri 1.80 `RR' (6) "Ratio = `RR' 95% CI `ci'", ms(i) mlabcolor(black))"'
			
		}
		
		if "`key'" == "" {
			local draw  "nodraw"
		}

		#delimit;
		tw 
			(area `delta_ly' `x', horizontal fcolor("`color1'") lcolor(white) fintensity(15)  base(`limitl')) 
			(area `delta_uy' `x', horizontal fcolor("`color2'") lcolor(white) fintensity(15)   base(`delta_l')) 
			(area `limituy' `x', horizontal fcolor("`color3'") lcolor(white)  fintensity(15)  base(`delta_u'))
			(scatteri 2 `newl' (0) "`syml'", ms(i) mlabcolor(black))
			(scatteri 2 `newu' (0) "`symu'", ms(i) mlabcolor(black))
			(scatteri 3 `midp1' (6) "Inferiority", ms(i) mlabcolor(black))
			(scatteri 3 `midp2' (6) "Similarity", ms(i) mlabcolor(black))
			(scatteri 3 `midp3' (6) "Superiority", ms(i) mlabcolor(black))
			(pci 2 `newl' 2 `newu', lcolor(black))
			`showstat'
			`vline'
			,
			`options'
			`draw'
			name(ci, replace)
		;
		#delimit cr
		
		if "`key'" == "" {
			qui {
				clear
				set obs 6
				egen y = seq(), f(1) t(6)

				forvalues r=1/6{
					gen y`r' = y[`r']/2
				}
				replace y = y*0.25
				gen y7 = y + 0.25

				drop y
				egen y = fill(50 50 100 100 150 150)
				replace y = y*0.01
			}
			
			local forbidden "fysize nodraw name"
			foreach opt of local forbidden {
				if strpos(`"`keyoptions'"', substr(`"`opt'"', 1, 4)) != 0 {
					di as error "Option `opt'() not allowed"
					exit
				}
			}
		
			local myopts `""yscale(noextend noline)""' 
			local myopts `"`myopts' `"ylabel("", noticks nolabels gextend nogrid)"'"'
			local myopts `"`myopts' `"legend(off)"'"'
			local myopts `"`myopts' `"graphregion(margin(r=0 l=0) color(white))"'"'
			local myopts `"`myopts' `"xlabel("", noticks nolabels gextend axis(1))"'"'
			local myopts `"`myopts' `"xscale(noline axis(1))"'"'
			local myopts `"`myopts' `"title("Key:", place(9))"'"'
			local myopts `"`myopts' `"xtitle(" ")"'"'
			local myopts `"`myopts' `"subtitle("    Condition					  			Conclusion", place(9))"'"'
			
			foreach opt of local myopts {
				if strpos(`"`keyoptions'"', substr(`"`opt'"', 1, 4)) == 0 {
					local keyoptions `"`keyoptions' `opt'"'
				}
			}
				
			#delimit;
			tw 	
				(area y2 y, fcolor("gray") lcolor(white) fintensity(15)  base(0.5))
				(area y3 y, fcolor("red") lcolor(white) fintensity(15)  base(1))
				(area y4 y, fcolor("orange") lcolor(white) fintensity(15)  base(1.5))
				(area y5 y, fcolor("green") lcolor(white) fintensity(15)  base(2))
				(area y6 y, fcolor("gray") lcolor(white) fintensity(15)  base(2.5))
				(scatter y7 y1, ms(i) 
					text(2.75 0.5 "Upper CI < `delta_uorig'  &  Lower CI < `delta_l'", place(e))
					text(2.25 0.5 "Lower CI >= `delta_uorig'", place(e))
					text(1.75 0.5 "Upper CI <= `delta_uorig'  &  Lower CI >= `delta_l'", place(e))
					text(1.25 0.5 "Upper CI < `delta_l'", place(e))
					text(0.75 0.5 "Upper CI > `delta_uorig'  &  Lower CI > `delta_l'", place(e))
					
					text(2.75 1.25 "Non-superiority", place(e))
					text(2.25 1.25 "Superiority", place(e))
					text(1.75 1.25 "Similarity", place(e))
					text(1.25 1.25 "Inferiority", place(e))
					text(0.75 1.25 "Non-inferiority", place(e)))
				(pci 3 1.1 0.5 1.1, lcolor(white) lw(vvthick))
					,
					`keyoptions'
					fysize(120)
					nodraw
					name(key, replace)
			;
			#delimit cr
			
			#delimit;
			tw 
				(scatteri 1 1, ms(i) mlabcolor(white))
				(scatteri 2 2, ms(i) mlabcolor(white))
				,
				graphregion(margin(r=0 l=0) color(white))
				plotregion(color(white))
				yscale(noextend noline)
				xscale(noline axis(1))
				legend(off)
				xlabel("", labc(none) tlcolor(none) gextend axis(1))
				ylabel(, noticks nolabels gextend axis(1) nogrid)
				xtitle(" ")
				ytitle(" ")
				fysize(30)
				nodraw
				name(blankg, replace)
			;
			#delimit cr	
			
			#delimit ;
				graph combine blankg key blankg
				, 
				cols(1)
				graphregion(margin(r=0 l=0) color(white))
				name(key, replace)
				imargin(0 0 0)
				;
			#delimit cr


			#delimit ;
				graph combine ci key
				, 
				rows(1)
				graphregion(margin(r=0 l=0) color(white))
				name(key, replace)
				imargin(0 0)
				;
			#delimit cr
		}
		restore
	}
	
end


/*
cap program drop ci3
program define ci3
version 14.0

end
*/

cap program drop zvalue
program define zvalue, rclass
version 14.0

	syntax anything(name=data id="data"), delta(string)
	
	tokenize `data'
	local p1 = `1'
	local p0 = `2'
	local p10 = `3'
	local p01 = `4'
	local p00 = `5'
	local n = `6'
	
	local p10_tilde = (-`p1' + `delta'^2 *(`p0' + 2*`p10') + sqrt((`p1' - `delta'^2*`p0')^2 + 4* `delta'^2 * `p10' * `p01') )/(2*`delta'*(`delta' + 1))
    local p01_tilde = `delta' * `p10_tilde' - (`delta' - 1)*(1 - `p00')
	
	local z = (sqrt(`n')*(`p1' - `delta'*`p0'))/sqrt(`delta'*(`p10_tilde' + `p01_tilde'))
	
	return local z `z'
end


/*Iterative procedure*/
version 14.0
mata:
	mata clear

	void cml_ci(real rowvector x, real scalar alpha) {
	
	real scalar a, b, c, d, n, p1, p0, RR, inits_l, inits_u

	ci = J(2, 2, .)

	a = x[1]
	b = x[2]
	c = x[3]
	d = x[4]
	
	n = a + b + c + d
	
	p1 = (a + b)/n
	p0 = (a + c)/n
	
	RR = p1/p0
	
	level = (alpha*2, alpha)
	
	for (i=1; i<=2; i++) {
	
	inits_l = RR*exp(-1*invnormal(1 - level[i]/2)*sqrt((b + c)/((a + b)*(a + c))))
	inits_u = RR*exp(invnormal(1 - level[i]/2)*sqrt((b + c)/((a + b)*(a + c))))	
		
	S = solvenl_init()
	solvenl_init_evaluator(S, &fun())
	solvenl_init_type(S, "zero")
	solvenl_init_technique(S, "newton")
	solvenl_init_numeq(S, 1)
	solvenl_init_narguments(S, 1)	
	
	solvenl_init_argument(S, 1, (a, b, c, n, alpha*2/i, 1))
	solvenl_init_startingvals(S, inits_l)
	solvenl_init_iter_log(S, "off")
	lower = solvenl_solve(S)

	solvenl_init_argument(S, 1, (a, b, c, n, alpha*2/i, -1))
	solvenl_init_startingvals(S, inits_u)
	solvenl_init_iter_log(S, "off")
	upper = solvenl_solve(S)
	
	ci[i, .] = (lower, upper)
	}
	st_matrix("ci", ci)
}
end

version 14.0
mata:	
	void fun(real scalar RR, real scalar y, real rowvector x){
		a = x[1]
		b = x[2]
		c = x[3]
		n = x[4]
		level = x[5]
		sign = x[6]
		
		real scalar p10_tilde
		real scalar p01_tilde
		
		z_alpha = invnormal(level/2)
		
		A = n*(1 + RR)
		B = (c + a)*RR^2 - (a + b + 2*c)
		C = c*(1 - RR)*(a + b + c)/n
		p01_tilde = (-B + sqrt(B^2 - 4*A*C))/(2*A)
		p10_tilde = ((a + b + c)/n - p01_tilde)/RR
		y = ((b + a - (a + c)*RR)/(RR*(2*p01_tilde + p10_tilde*(RR - 1))))/sqrt((n)/(RR*(2*p01_tilde + p10_tilde*(RR - 1)))) + sign*z_alpha
	}
end
