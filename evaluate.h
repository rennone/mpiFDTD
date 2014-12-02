#ifndef _EVALATE_H__
#define _EVALATE_H__

//各波長の角反射率
// reflec[lambda][degree] : 波長λのdegree°方向の反射率
extern double evaluate_evaluate(double **reflec, int stLambda, int enLambda);

#endif
