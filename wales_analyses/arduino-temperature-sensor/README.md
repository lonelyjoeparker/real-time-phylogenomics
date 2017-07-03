# Arduino temperature monitoring script

## Motivation

A simple temperature sensor used to monitor the internal conditions of polyboxes used for field transportation and storage of temperature-sensitive reagents in Wales.

## Credit

This code is modified (a bit) from learningaboutelectronics.com’s [example](http://www.learningaboutelectronics.com/Articles/LM335-temperature-sensor-circuit.php)

## Requirements

 - Arduino Uno / Genuino
 - Wire (~5mm dia gives good conductance for electricity without seeming to cold-bridge to the outside of the box ‘too’ badly - based on ‘feeling the wire next to a -20ºC box’ tests), solder, tape etc
 - (LM335AZ)[https://www.google.co.uk/search?q=lm335az&oq=lm335az&aqs=chrome..69i57.2262j0j7&sourceid=chrome&ie=UTF-8] temperature transducer and 2KOhm resistor.
 - 4x standard LEDs and 1kOhm resistors.

## Setup / construction

Pinouts as for the example code with LEDs in pins 4 (red), 7 & 8 (orange) and green (12). Additionally the on-board (Uno) pin 12 LED is used for the blink function (visually verifies the script hasn’t hung.)

To set up the transducers, solder the connections to one end of the wire as in the examples and product specification, leaving the other wire end blunt. Use an awl or similar tool to pierce a hole in the polybox of choice 50mm<depth<150mm from the top of the box. You want to make the smallest possible hole to keep a decent heat seal after inserting the wire. Push the wire through, prune and connect to the Arduino checking it works (a lab thermometer and some icy water are handy here) then solder/croc clip as desired (we used one Arduino to read several boxes, powered via USB or 9v).

## Usage

The LEDs in combination will tell you if the box temperature is:
 - below ~10ºC (green LED)
 - above ~6ºC (red LED)
 - decreased since last measurement (orange LED 1)
 - increased since last measurement (orange LED 2).
 - whether the script is still responsive (on-board orange LED).

In combination these will let you work out whether box temperature is suitable for prepared libraries and flowcells (no red, no green); for storage DNA and kits/enzymes (green, no red) or too warm (red).