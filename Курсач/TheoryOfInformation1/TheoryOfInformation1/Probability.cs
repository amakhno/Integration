using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TheoryOfInformation1
{
    class Probability
    {
        double[] ArrayOfProbabilityA;
        double[] ArrayOfProbabilityB;
        double[,] Matrix;

        public double AverageTogetherInformation;

        double Ha, Hb, H_ab, H_ba;
        double H_AB;//Энтропии H(x) H(y) H(X|Y) H(Y|X) H(XY)

        int Size = 2;
        int SizeN;
        int SizeM;

        public Probability()
        {
            ArrayOfProbabilityA = new double[] { 1/3.0, 1/3.0 , 1/3.0};
            ArrayOfProbabilityB = new double[] { 2/3.0, 1/3.0 };
            /*
            ArrayOfProbabilityA = new double[] { 0.25, 0.69, 0.06 };
            ArrayOfProbabilityB = new double[] { 0.51, 0.48, 0.01 };

            
            Random rand = new Random();
            Matrix = new double[Size, Size];
            for (int i = 0; i < Size; i++)
                for (int j = 0; j < Size; j++)
                    Matrix[i, j] = ArrayOfProbabilityB[i] * ArrayOfProbabilityA[j];
            */
            Matrix = new double[ArrayOfProbabilityB.Count(), ArrayOfProbabilityA.Count()];
            /*Matrix[0, 0] = 1/3.0;
            Matrix[0, 1] = 1 / 3.0;
            Matrix[1, 0] = 0.0;
            Matrix[1, 1] = 1/3.0;*/
            Matrix[0, 0] = 0;
            Matrix[0, 1] = 1 / 3.0;
            Matrix[0, 2] = 1 / 3.0;
            Matrix[1, 0] = 1 / 3.0;
            Matrix[1, 1] = 0;
            Matrix[1, 2] = 0;
            SizeN = ArrayOfProbabilityA.Count();
            SizeM = ArrayOfProbabilityB.Count();
            AverageTogetherInformation = CalculateAverageTogetherInformation();
            EntropyCalc();
        }

        public double GetEntropyOfA()
        {
            double Ha = 0;
            for (int i = 0; i < ArrayOfProbabilityA.Count(); i++)
                Ha -= ArrayOfProbabilityA[i] * Math.Log(ArrayOfProbabilityA[i], 2);
            this.Ha = Ha;
            return Ha;
        }

        public double GetEntropyOfB()
        {
            double Ha = 0; ;
            for (int i = 0; i < ArrayOfProbabilityB.Count(); i++)
                Ha -= ArrayOfProbabilityB[i] * Math.Log(ArrayOfProbabilityB[i], 2);
            this.Hb = Ha;

            return Ha;
        }

        public double GetEntropyAIfB()
        {
            return H_ab;
        }

        public double GetEntropyBIfA()
        {
            return H_ba;
        }

        public double GetEntropyAWithB()
        {
            return H_AB;
        }

        public void EntropyCalc()
        {
            Ha = GetEntropyOfA();

            Hb = GetEntropyOfB();

            H_ab = Ha - AverageTogetherInformation;

            H_ba = Hb - AverageTogetherInformation;

            H_AB = Ha + Hb - AverageTogetherInformation;
        }

        private double CalculateAverageTogetherInformation()
        {
            double AverMutInf;

            double[,] MutInfMat = new double[SizeM, SizeN]; 

            double sumA, sumB;

            for (int i = 0; i < SizeM; i++)
            {
                for (int j = 0; j < SizeN; j++)
                {
                    sumA = 0; sumB = 0;
                    for (int k = 0; k < SizeM; k++)
                        sumA += Matrix[k, j];

                    for (int l = 0; l < SizeN; l++)
                        sumB += Matrix[i, l];

                    if (Matrix[i, j] != 0)
                        MutInfMat[i, j] = Math.Log(Matrix[i, j] / (sumA * sumB), 2);
                    else
                        MutInfMat[i, j] = 0;
                }
            }

            AverMutInf = 0;
            for (int i = 0; i < SizeM; i++)
                for (int j = 0; j < SizeN; j++)
                    AverMutInf += Matrix[i, j] * MutInfMat[i, j];

            return AverMutInf;
        }

        public double GetAverageTogetherInformation()
        {
            return AverageTogetherInformation;
        }
    }
}
