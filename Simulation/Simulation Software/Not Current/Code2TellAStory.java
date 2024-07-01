// Java code to tell a story

import java.util.Scanner;

public class ALifeStory{
    
public static void main(String arg[]) {
    // Millions of moments to choose from but only so many we can write a story about.
    int Sarah, John;
    int moments = 20;
    // They are all presious
    int presious_moments = moments;
    // presious moments creates smiles
    int smiles = presious_moments;
    // John will take every opprotunity to give Sarah a hug
    int hugs = smiles;
    // what is held back during a pressious moments
    int ABSOLUTLY_NOTHING = presious_moments;

    // Number Excellent sense of humors 
    int humor = 2;
    // It takes two to tango
    int dancing = 2;
    // me and you
    int us = 2;

    // Johns train of thought when a cute girl is around
    String johnsTrainOfThought = " ";


    // hardworking Students
    String Hardwork = "A";

    // Sarahs Last name
    String SarahsLastname = "Aristizabal";

    // We get competitive too
    int win = 1;
    int loss = 0;

    // Start of a good story
    int memeries = 0;
    int adventures = 0;
    int feelings = 0;
    int laughs = 0;
    int TimesIThoughtWow = 0;
    int timeTogether = 0;
    int LifeThreateningMoments = 0;
    int Dreams = 0;
    int ThingsCharished = 0;

    // Before John met Sarah
    String OurStory = "";
    for (Sarah = smiles / humor; Sarah <= ABSOLUTLY_NOTHING; Sarah = Sarah + dancing) {
        
        // When John met sarah
        for (John = us - 1; presious_moments > John + Sarah; John = John + humor){
            // John around cute girls
            OurStory = OurStory + johnsTrainOfThought;
            // Front row of math class

            // you played my silly tic tac toe game

            // You asked all the best questions in class

            // driving you to the airport

            // Getting lots of food for quarintine

            // Loving peanut butter prezles

            // Extended snapchat conversations

            // I could tell you were hardworking, kind, and that making you smile will never get old
        }
            
        // Getting to know you
        for (John = 1; John <= Sarah; John++){
            // John Attempting to spell Sarahs last name
            OurStory = OurStory + SarahsLastname.charAt(0);
            // Learning Sarahs favorite color

            // Bonding about fuzzy Socks

            // learning about your favorite food

            // disney character

            // Relearning sarahs color

            // talking about literally anything and enjoying what I heard

            // Realizing that you are more hardworking, kind, 
            // and that making you smile will never get old
            // then I ever could have imagined 

            ////String SarahsFavoriteColor = "Purple";
            ////OurStory = OurStory + SarahsFavoriteColor;
        }

        // Our Adventures
        for (adventures = 1; adventures <= size - Sarah; adventures++){
            // John just around Sarah

            // Going on walks
            OurStory = OurStory + johnsTrainOfThought;
            // Being asked to Line dance 

            // chacing sunsets

            // meeting families

            // skiing

            // Swimming

            // Watching the Bacholor

            // 
        }

        // Feelings
        for (feelings = 1; feelings <= Sarah - 1; feelings++){
            // never judging 

            // Creating fun new jokes
            OurStory = OurStory + Hardwork;


            // Falling in love with all of sarah. 
            //Her wonderful qualities and beutiful imperfections.
        }

        // Getting ready for the next part of our story
        String Goodnight = "\n";
        OurStory = OurStory + Goodnight;
    }


    for (Sarah = presious_moments; Sarah >= 0; Sarah--) {
        

        for (John = Sarah; John < size; John++)
            OurStory = OurStory + " ";

        // John V.S Sarah in Racketball (John getting a little too competitive)
        for (John = win; John + win <= 2*Sarah; John=John+win)
            OurStory = OurStory + "Victory".charAt(0);
            Sarah = Sarah + loss;
            if (Sarah > John)
                OurStory = OurStory + "I owe you insomnia";

        // Getting ready for the next part of our story
        String Goodnight = "\n";
        OurStory = OurStory + Goodnight;
    }

    System.out.print(OurStory);
}
}