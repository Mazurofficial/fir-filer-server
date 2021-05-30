import freqz from 'ndarray-freqz';
import zeros from "zeros";
import cops from 'ndarray-complex'

const calcKImpulseResponse = (params) =>{
    let Fs = params.Fs
    let Fa = params.Fa
    let Fb = params.Fb
    let o = params.order || 51
    let alpha = params.Att || 100
    const ino = (val) => {
        let d = 0
        let ds = 1
        let s = 1
        while (ds > s * 1e-6) {
            d += 2
            ds *= val * val / (d * d)
            s += ds
        }
        return s
    }
    if (o / 2 - Math.floor(o / 2) === 0) {
        o++
    }
    let Np = (o - 1) / 2
    let A = []
    let beta = 0
    let cnt = 0
    let inoBeta
    let ret = []

    A[0] = 2 * (Fb - Fa) / Fs
    for (cnt = 1; cnt <= Np; cnt++) {
        A[cnt] = (Math.sin(2 * cnt * Math.PI * Fb / Fs) - Math.sin(2 * cnt * Math.PI * Fa / Fs)) / (cnt * Math.PI)
    }
    // empirical coefficients
    if (alpha < 21) {
        beta = 0
    } else if (alpha > 50) {
        beta = 0.1102 * (alpha - 8.7)
    } else {
        beta = 0.5842 * Math.pow((alpha - 21), 0.4) + 0.07886 * (alpha - 21)
    }

    inoBeta = ino(beta)
    for (cnt = 0; cnt <= Np; cnt++) {
        ret[Np + cnt] = A[cnt] * ino(beta * Math.sqrt(1 - (cnt * cnt / (Np * Np)))) / inoBeta
    }
    for (cnt = 0; cnt < Np; cnt++) {
        ret[cnt] = ret[o - 1 - cnt]
    }
    return ret
}
const calcImpulseResponse = (params) => {
    let Fs = params.Fs
    let Fc = params.Fc
    let o = params.order
    let omega = 2 * Math.PI * Fc / Fs
    let cnt = 0
    let dc = 0
    let ret = []
    // sinc function is considered to be
    // the ideal impulse response
    // do an idft and use Hamming window afterwards
    for (cnt = 0; cnt <= o; cnt++) {
        if (cnt - o / 2 === 0) {
            ret[cnt] = omega
        } else {
            ret[cnt] = Math.sin(omega * (cnt - o / 2)) / (cnt - o / 2)
            // Hamming window
            ret[cnt] *= (0.54 - 0.46 * Math.cos(2 * Math.PI * cnt / o))
        }
        dc = dc + ret[cnt]
    }
    // normalize
    for (cnt = 0; cnt <= o; cnt++) {
        ret[cnt] /= dc
    }
    return ret
}
const invert = (h) => {
    let cnt
    for (cnt = 0; cnt < h.length; cnt++) {
        h[cnt] = -h[cnt]
    }
    h[(h.length - 1) / 2]++
    return h
}
const bs = (params) => {
    let lp = calcImpulseResponse({
        order: params.order,
        Fs: params.Fs,
        Fc: params.F2
    })
    let hp = invert(calcImpulseResponse({
        order: params.order,
        Fs: params.Fs,
        Fc: params.F1
    }))
    let out = []
    for (let i = 0; i < lp.length; i++) {
        out.push(lp[i] + hp[i])
    }
    return out
}


const fir = (filterType,params) => {
    let coefficients = [];
    switch (filterType){
        case 'LPF':
            coefficients = calcImpulseResponse(params);
            break;
        case 'HPF':
            coefficients = invert(calcImpulseResponse(params));
            break;
        case 'BSF':
            coefficients = bs(params);
            break;
        case 'BPF':
            coefficients = invert(bs(params));
            break;
        case 'KBF':
            coefficients = calcKImpulseResponse(params);
            break;
        default :
            return('no such method')
    }
    const freq = freqz(coefficients,[1])
    const f_re = freq.H_r.data
    const f_im = freq.H_i.data
    let result = []
    const magnitude = zeros([512])
    cops.mag(magnitude, freq.H_r, freq.H_i);
    const frequency = magnitude.data
    for (let i = 0; i <= params.Fs/2; i++){
        // result.push([0.5*params.Fs*f_re[i]/Math.PI,20*Math.log10(Math.abs(f_im[i]))])
        result.push([i, f_im[i*2]])
    }
    return result
}

export default fir

// console.log(fir('LPF',{
//     Fs: 1000,
//     order: 1000,
//     Fc: 100
// }))