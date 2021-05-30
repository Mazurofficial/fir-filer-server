import BESSEL from 'bessel';
import freqz from 'ndarray-freqz';

const isEven = n => {
    return n % 2 === 0;
}

const calcFcLPF = (Fa,Fp) => {
    const Bt = Fa - Fp;
    const Fc = Fp + Bt/2;
    return {Bt: Bt,Fc: Fc}
}
const calcFcHPF = (Fa,Fp) => {
    const Bt = Fp - Fa;
    const Fc = Fa + Bt/2;
    return {Bt: Bt,Fc: Fc}
}
const calcFcBPF = (Fa1,Fp1,Fa2,Fp2) => {
    const Bt = Math.min((Fp1 - Fa1),(Fa2-Fp2))
    const Fc1 = Fp1 - Bt/2;
    const Fc2 = Fp2 + Bt/2;
    return {Bt: Bt,Fc1: Fc1,Fc2: Fc2}
}
const calcFcBSF = (Fa1,Fp1,Fa2,Fp2) => {
    const Bt = Math.min((Fa1 - Fp1),(Fp2-Fa2))
    const Fc1 = Fp1 + Bt/2;
    const Fc2 = Fp2 - Bt/2;
    return {Bt: Bt,Fc1: Fc1,Fc2: Fc2}
}

const calcD = A =>{
    if (A <= 21){
        return 0.9222
    }else if (A > 21){
        return (A - 7.95)/14.36
    }
}
const calcAlpha = A =>{
    if (A <= 21){
        return  0
    } else if(A <= 50){
        return  0.5842*Math.pow(A - 21, 0.4) + 0.07886*(A -21)
    } else if(A > 50){
        return  0.1102*(A - 8.7)
    }
}
const calcFilterOrder = (D,Fs,Bt) => Math.round((Fs*D)/Bt)

const calcFourierCoefficientsNEven = (type, Fs, Fc, Fc1, Fc2, N) => {
    let a = [];
    switch (type){
        case 'LPF':
            for (let i = 0; i <= (N/2); i++){
                //a.push((1/(Math.PI*(i - 1/2))) * Math.sin(2*Math.PI*(i - 1/2)*(Fc/Fs)))
                a.push(Math.sin(2*Math.PI*(i - 0.5)*Fc/Fs)/(Math.PI*(i - 0.5)))
            }
            a.shift()
            break
        case 'HPF':
            for (let i = 0; i <= N/2; i++){
                a.push((1/(Math.PI*(i - 1/2))) * (Math.sin(Math.PI*(i - 1/2)) - Math.sin(2*Math.PI*(i - 1/2)*(Fc/Fs))))
            }
            a.shift()
            break
        case 'BPF':
            for (let i = 0; i <= N/2; i++){
                a.push((1/(Math.PI*(i - 1/2))) * (Math.sin(2*Math.PI*(i - 1/2)*(Fc2/Fs)) - Math.sin(2*Math.PI*(i - 1/2)*(Fc1/Fs))))
            }
            a.shift()
            break
        case 'BSF':
            for (let i = 0; i <= N/2; i++){
                a.push((1/(Math.PI*(i - 1/2))) * (Math.sin(Math.PI*(i - 1/2)) - Math.sin(2*Math.PI*(i - 1/2)*(Fc2/Fs)) + Math.sin(2*Math.PI*(i - 1/2)*(Fc1/Fs))))
            }
            a.shift()
            break
        default:
            return 'No such method'
    }
    return a
}
const calcFourierCoefficientsNOdd = (type, Fs, Fc, Fc1, Fc2, N) => {
    let a = [];
    switch (type){
        case 'LPF':
            a.push(2*(Fc/Fs))
            for (let i = 1; i <= (N+1)/2; i++){
                // a.push((1/(Math.PI*i)) * Math.sin(2*Math.PI*i*(Fc/Fs)))
                a.push(Math.sin(2*Math.PI*i*Fc/Fs)/(Math.PI*i))
            }
            // a.shift()
            a.pop()
            break
        case 'HPF':
            a.push(1 - 2*(Fc/Fs))
            for (let i = 1; i <= N/2; i++){
                a.push(-(1/(Math.PI*i)) * Math.sin(2*Math.PI*i*(Fc/Fs)))
            }
            a.pop()
            break
        case 'BPF':
            a.push(2*((Fc2-Fc1)/Fs))
            for (let i = 1; i <= N/2; i++){
                a.push((1/(Math.PI*i)) * (Math.sin(2*Math.PI*i*(Fc2/Fs)) - Math.sin(2*Math.PI*i*(Fc1/Fs))))
            }
            a.pop()
            break
        case 'BSF':
            a.push(1 - 2*((Fc2-Fc1)/Fs))
            for (let i = 1; i <= N/2; i++){
                a.push((1/(Math.PI*i)) * (Math.sin(2*Math.PI*i*(Fc1/Fs)) - Math.sin(2*Math.PI*i*(Fc2/Fs))))
            }
            a.pop()
            break
        default:
            return 'No such method'
    }
    return a
}

const calcKaizerWindow = (alpha, M) => {
    let w = []
    for (let i = 0; i <= M/2; i++){
        let beta = alpha*Math.sqrt(1-Math.pow((2*i)/M,2))
        w.push(BESSEL.besseli(beta,0)/BESSEL.besseli(alpha,0))
    }
    return w
}
const calcResultCoefficients = (a, w, M) => {
    let resultArr = [];
    for (let i = 0; i <= M/2; i++){
        resultArr.push(a[i] * w[i])
    }
    return resultArr
}
const calcFirFilterCoefficients = (coeffs, N) => {
    let h = [...coeffs]
    let hNew = [...h]
    let hMirror = []
    if(isEven(N)){
        hMirror = h.reverse()
    }else {
        h.reverse().pop()
        hMirror = h
    }
    return hMirror.concat(hNew)

}

const calcFR = (N,FS, coeffs) => {

    const freq = freqz(coeffs,[1])
    const f_re = freq.H_r.data
    const f_im = freq.H_i.data
    let frequencyResponse = []
    for (let i = 0; i < FS/2; i++){
        frequencyResponse.push(Math.sqrt(Math.pow(f_re[i],2) + Math.pow(f_im[i],2)))
    }
    let f =[]
    for (let i = 0; i < 500; i++ ){
        f.push(i)
    }
    let result = []
    for (let i = 0; i <= f.length; i++){
        //result.push([i,20*Math.log10(frequencyResponse[i])])
        result.push([i,frequencyResponse[i]])
    }
    return result
}

const calcCoeffsGeneral = (type,o,Fs,Bt,Fc,Fc1,Fc2,N) => {
    const A = -20*Math.log10(o);
    const D = calcD(A);
    const alpha = calcAlpha(A);
    let Order
    let M
    if (N !== 0){
        Order = N
        M = Order -1;
    }else {
        M = calcFilterOrder(D,Fs,Bt);
        Order  = M + 1;
    }
    let a =[];
    if (isEven(Order)){
        a = calcFourierCoefficientsNEven(type,Fs,Fc,Fc1,Fc2,Order)
        console.log('even')
    }else {
        a = calcFourierCoefficientsNOdd(type,Fs,Fc,Fc1,Fc2,Order)
        console.log('noteven')
    }
    const w = calcKaizerWindow(alpha,M)
    let coefficients = calcResultCoefficients(a,w,M)
    console.log('a')
    console.log(a)
    const coeffs = calcFirFilterCoefficients(a, Order)
    const kaizerCoeffs = calcFirFilterCoefficients(coefficients, Order)
    console.log('Standart coeffs')
    console.log(coeffs)
    console.log('Kaizer coeffs')
    console.log(kaizerCoeffs)
    return {
        coeffsFR: calcFR(Order, Fs, coeffs),
        kaizerCoeffsFR: calcFR(Order, Fs, kaizerCoeffs)
    }
}

const lpf = {
    type: 'LPF',
    Fc: 100,
    N: 10,
    Ap: 2,
    Aa: 30,
    Fs: 1000
}

const hpf = {
    type: 'HPF',
    Fp: 150,
    Fa: 100,
    Ap: 2,
    Aa: 30,
    Fs: 1000
}

const bpf = {
    type: 'BPF',
    Fc2: 400,
    Fc1: 300,
    Ap: 2,
    Aa: 30,
    Fs: 1000,
    N: 10
}

const bsf = {
    type: 'BSF',
    Fp1: 130,
    Fa1: 220,
    Fp2: 250,
    Fa2: 170,
    Ap: 2,
    Aa: 30,
    Fs: 1000
}

const firFilter = (data) => {
    let result;
    const o1 = (Math.pow(10, 0.05 * data.Ap) + 1) / (Math.pow(10, 0.05 * data.Ap) - 1);
    const o2 = Math.pow(10, -0.05 * data.Aa);
    const o = Math.min(o1,o2);
    console.log(data)
    switch (data.type){
        case 'LPF':
            if (data.Fc === undefined){
                const fc_bt = calcFcLPF(data.Fa,data.Fp)
                result = calcCoeffsGeneral(data.type,o,data.Fs,fc_bt.Bt,fc_bt.Fc,0,0,0)
            }else {
                console.log('i recieve data')
                result = calcCoeffsGeneral(data.type,o,data.Fs,0,data.Fc,0,0,data.N)
            }
            break
        case 'HPF':
            if (data.Fc === undefined){
                const fc_bt = calcFcHPF(data.Fa,data.Fp)
                result = calcCoeffsGeneral(data.type,o,data.Fs,fc_bt.Bt,fc_bt.Fc,0,0,0)
            }else {
                result = calcCoeffsGeneral(data.type,o,data.Fs,0,data.Fc,0,0,data.N)
            }
            break
        case 'BPF':
            if (data.Fc1 === undefined){
                const fc_bt = calcFcBPF(data.Fa1,data.Fp1,data.Fa2,data.Fp2)
                result = calcCoeffsGeneral(data.type,o,data.Fs,fc_bt.Bt,0,fc_bt.Fc1,fc_bt.Fc2,0)
            }else {
                result = calcCoeffsGeneral(data.type,o,data.Fs,0,0,data.Fc1,data.Fc2,data.N)
            }
            break
        case 'BSF':
            if (data.Fc1 === undefined){
                const fc_bt = calcFcBSF(data.Fa1,data.Fp1,data.Fa2,data.Fp2)
                result = calcCoeffsGeneral(data.type,o,data.Fs,fc_bt.Bt,0,fc_bt.Fc1,fc_bt.Fc2,0)
            }else {
                result = calcCoeffsGeneral(data.type,o,data.Fs,0,0,data.Fc1,data.Fc2,data.N)
            }
            break
        default:
            return 'NO SUCH METHOD'
    }
    return result
}

export default firFilter
console.log(firFilter(bpf))
