using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace TheoryOfInformation1
{
    class Program
    {
        static void Main(string[] args)
        {
            Probability prob = new Probability();

            Console.WriteLine("Энтропия ансамбля А: " + prob.GetEntropyOfA().ToString());
            Console.WriteLine("Энтропия ансамбля В: " + prob.GetEntropyOfB().ToString());

            Console.WriteLine("Средняя взаимная информация " + prob.GetAverageTogetherInformation().ToString());

            Console.WriteLine("Энтропия A|B " + prob.GetEntropyAIfB().ToString());
            Console.WriteLine("Энтропия B|A " + prob.GetEntropyBIfA().ToString());

            Console.WriteLine("Совместная энтропия " + prob.GetEntropyAWithB().ToString());

            Console.ReadKey();
        }
    }
}
