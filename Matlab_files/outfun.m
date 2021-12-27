function stop=outfun(x,optimValues,state)
 load('history')
   stop=false;
           history.fval = [history.fval; optimValues.fval];
           history.x = [history.x x];
           history.gradient=[history.gradient optimValues.gradient];
           save('history','history');

end
