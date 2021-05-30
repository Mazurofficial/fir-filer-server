import fir from "./fir.js";
import firFilter from './firfilter.js'
import express from "express";
import cors from "cors";

const app = express();
const router = express.Router();
app.use(cors());
app.use(express.json());
app.use("/", router);
app.listen(process.env.PORT || 5000, () => console.log("Server Running"));

router.post("/", (req, res) => {
    let firFilterCoeffs;
    const filterType = req.body.filterType;
    const Fp = req.body.Fp;
    const Fa = req.body.Fa;
    const Ap = req.body.Ap;
    const Aa = req.body.Aa;
    const Fs = req.body.Fs;
    const Fp1 = req.body.Fp1;
    const Fp2 = req.body.Fp2;
    const Fa1 = req.body.Fa1;
    const Fa2 = req.body.Fa2;
    const Fc = req.body.Fc;
    const Fc1 = req.body.Fc1;
    const Fc2 = req.body.Fc2;
    const N = req.body.N;
    if(N !== undefined){
        if (Fc !== undefined){
            firFilterCoeffs = firFilter(
                {
                    type: filterType,
                    Fc: Fc,
                    Ap: Ap,
                    Aa: Aa,
                    Fs: Fs,
                    N: N
                });
        } else {
            firFilterCoeffs = firFilter(
                {
                    type: filterType,
                    Fc1: Fc1,
                    Fc2: Fc2,
                    Ap: Ap,
                    Aa: Aa,
                    Fs: Fs,
                    N: N
                });
        }
    }else{
        if (Fp !== undefined){
            firFilterCoeffs = firFilter(
                {
                    type: filterType,
                    Fp: Fp,
                    Fa: Fa,
                    Ap: Ap,
                    Aa: Aa,
                    Fs: Fs
                });
        } else {
            firFilterCoeffs = firFilter(
                {
                    type: filterType,
                    Fp1: Fp1,
                    Fa1: Fa1,
                    Fp2: Fp2,
                    Fa2: Fa2,
                    Ap: Ap,
                    Aa: Aa,
                    Fs: Fs
                });
        }
    }
    let coeffsFR = Array.from(firFilterCoeffs.coeffsFR)
    let kaizerCoeffsFR = Array.from(firFilterCoeffs.kaizerCoeffsFR)
    // console.log(coeffsFR)
    // console.log(kaizerCoeffsFR)
    res.json({ coeffsFR: coeffsFR, kaizerCoeffsFR: kaizerCoeffsFR })
});