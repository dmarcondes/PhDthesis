///////////////////////////////////////////////////////////////////////////////////////////////
//////////A data-driven systematic, consistent and feasible approach to Model Selection////////
//////////4.3 Multilayer W-operator                                                    ////////
//////////PhD Thesis, Diego Marcondes                                                  ////////
//////////Universisty of São Paulo, 2022                                               ////////
///////////////////////////////////////////////////////////////////////////////////////////////

#ifndef PROGRESS_H_INCLUDED
#define PROGRESS_H_INCLUDED

using System;

void ProgressBar(int progress, int total)
    {

        //draw empty progress bar
        Console.CursorLeft = 0;
        Console.Write("["); //start
        Console.CursorLeft = 32;
        Console.Write("]"); //end
        Console.CursorLeft = 1;
        float onechunk = 30.0f / total;

        //draw filled part
        int position = 1;
        for (int i = 0; i < onechunk * progress; i++)
        {
            Console.BackgroundColor = ConsoleColor.Green;
            Console.CursorLeft = position++;
            Console.Write(" ");
        }

        //draw unfilled part
        for (int i = position; i <= 31; i++)
        {
            Console.BackgroundColor = ConsoleColor.Black;
            Console.CursorLeft = position++;
            Console.Write(" ");
        }

        //draw totals
        Console.CursorLeft = 35;
        Console.BackgroundColor = ConsoleColor.Black;
        Console.Write(progress.ToString() + " of " + total.ToString() + "    ");
    }

#endif // PROGRESS_H_INCLUDED
