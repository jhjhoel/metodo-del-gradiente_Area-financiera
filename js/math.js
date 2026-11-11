/* Núcleo numérico SIFIN Bolivia. Implementaciones propias para matrices densas pequeñas. */
const Numeric = (() => {
  const dot = (a,b) => a.reduce((s,v,i)=>s+v*b[i],0);
  const norm = a => Math.sqrt(dot(a,a));
  const matVec = (A,x) => A.map(row=>dot(row,x));
  const transpose = A => A[0].map((_,j)=>A.map(row=>row[j]));
  const identity = n => Array.from({length:n},(_,i)=>Array.from({length:n},(_,j)=>i===j?1:0));
  const subtract = (a,b) => a.map((v,i)=>v-b[i]);
  const scale = (a,s) => a.map(v=>v*s);
  const addScaled = (a,b,s) => a.map((v,i)=>v+s*b[i]);
  function cholesky(A){
    const n=A.length,L=Array.from({length:n},()=>Array(n).fill(0));
    for(let i=0;i<n;i++) for(let j=0;j<=i;j++){
      let sum=A[i][j]; for(let k=0;k<j;k++) sum-=L[i][k]*L[j][k];
      if(i===j){if(!(sum>1e-14)) throw new Error('La matriz no es definida positiva (Cholesky encontró un pivote no positivo).'); L[i][j]=Math.sqrt(sum);}
      else L[i][j]=sum/L[j][j];
    }
    return L;
  }
  function solveCholesky(A,b){
    const L=cholesky(A),n=b.length,y=Array(n).fill(0),x=Array(n).fill(0);
    for(let i=0;i<n;i++){let s=b[i];for(let j=0;j<i;j++)s-=L[i][j]*y[j];y[i]=s/L[i][i];}
    for(let i=n-1;i>=0;i--){let s=y[i];for(let j=i+1;j<n;j++)s-=L[j][i]*x[j];x[i]=s/L[i][i];}
    return x;
  }
  function jacobiEigenvalues(input){
    const A=input.map(r=>r.slice()),n=A.length,max=100*n*n;
    for(let it=0;it<max;it++){
      let p=0,q=1,maxVal=0;
      for(let i=0;i<n;i++)for(let j=i+1;j<n;j++)if(Math.abs(A[i][j])>maxVal){maxVal=Math.abs(A[i][j]);p=i;q=j;}
      if(maxVal<1e-12) break;
      const phi=.5*Math.atan2(2*A[p][q],A[q][q]-A[p][p]),c=Math.cos(phi),s=Math.sin(phi);
      const app=c*c*A[p][p]-2*s*c*A[p][q]+s*s*A[q][q];
      const aqq=s*s*A[p][p]+2*s*c*A[p][q]+c*c*A[q][q];
      for(let k=0;k<n;k++)if(k!==p&&k!==q){const akp=A[k][p],akq=A[k][q];A[k][p]=A[p][k]=c*akp-s*akq;A[k][q]=A[q][k]=s*akp+c*akq;}
      A[p][p]=app;A[q][q]=aqq;A[p][q]=A[q][p]=0;
    }
    return A.map((r,i)=>r[i]).sort((a,b)=>a-b);
  }
  function conditionNumber(A){const e=jacobiEigenvalues(A);if(e[0]<=0)return Infinity;return e[e.length-1]/e[0];}
  function isSymmetric(A,tol=1e-9){for(let i=0;i<A.length;i++)for(let j=i+1;j<A.length;j++)if(Math.abs(A[i][j]-A[j][i])>tol)return false;return true;}
  function validate(A,b){
    if(!Array.isArray(A)||A.length<4||A.some(r=>r.length!==A.length))throw new Error('La matriz debe ser cuadrada y tener al menos 4 filas y 4 columnas.');
    if(!Array.isArray(b)||b.length!==A.length)throw new Error('El vector b debe tener la misma dimensión que A.');
    if(A.flat().concat(b).some(v=>!Number.isFinite(v)))throw new Error('Todos los coeficientes deben ser números finitos.');
    if(!isSymmetric(A))throw new Error('Matriz no válida: A no es simétrica. Verifica aᵢⱼ = aⱼᵢ.');
    cholesky(A);
    if(norm(b)===0)throw new Error('El vector b no puede ser nulo para calcular el residuo relativo.');
    return true;
  }
  function conjugateGradient(A,b,x0,tol,maxIter,preconditioned=false){
    const n=b.length;let x=x0.slice(),r=subtract(b,matVec(A,x));
    const applyM = v => preconditioned?v.map((z,i)=>z/A[i][i]):v.slice();
    let z=applyM(r),p=z.slice(),rz=dot(r,z),rel=norm(r)/norm(b),iterations=0;
    const history=[rel];
    while(rel>=tol && iterations<maxIter){
      const Ap=matVec(A,p),den=dot(p,Ap);
      if(!(den>0)||!Number.isFinite(den))break;
      const alpha=rz/den;x=addScaled(x,p,alpha);r=addScaled(r,Ap,-alpha);
      rel=norm(r)/norm(b);iterations++;history.push(rel);
      if(rel<tol)break;
      z=applyM(r);const rzNew=dot(r,z);if(!Number.isFinite(rzNew)||Math.abs(rz)<1e-30)break;
      const beta=rzNew/rz;p=addScaled(z,p,beta);rz=rzNew;
    }
    const trueResidual=norm(subtract(b,matVec(A,x)))/norm(b);
    return {x,iterations,residual:trueResidual,converged:trueResidual<tol,history};
  }
  function preconditionedCondition(A){
    // M = diag(A). La matriz simétrica equivalente a M^-1 A es D^-1/2 A D^-1/2.
    const d=A.map((r,i)=>Math.sqrt(r[i]));
    const B=A.map((r,i)=>r.map((v,j)=>v/(d[i]*d[j])));
    return conditionNumber(B);
  }
  return {dot,norm,matVec,transpose,cholesky,solveCholesky,jacobiEigenvalues,conditionNumber,isSymmetric,validate,conjugateGradient,preconditionedCondition};
})();
