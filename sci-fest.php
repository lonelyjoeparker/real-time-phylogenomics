<!DOCTYPE html>
<html lang="en">
  <head>
	<!--
		UNCOMMENT TO REFRESH
    -->
	<meta http-equiv="refresh" content="8">
    <meta charset="utf-8">
    <meta http-equiv="X-UA-Compatible" content="IE=edge">
    <meta name="viewport" content="width=device-width, initial-scale=1">
    <!-- The above 3 meta tags *must* come first in the head; any other head content must come *after* these tags -->
    <meta name="description" content="">
    <meta name="author" content="">
    <!--<link rel="icon" href="../../favicon.ico">-->

    <title>Real-time phylogenomics</title>

    <!-- Bootstrap core CSS -->
    <link href="../css/bootstrap.css" rel="stylesheet">

	<!-- DataTables -->
	<!--<link rel="stylesheet" type="text/css" href="//cdn.datatables.net/1.10.11/css/jquery.dataTables.css">-->
	<link rel="stylesheet" type="text/css" href="../DataTables/DataTables-1.10.11/css/jquery.dataTables.css">

    <!-- HTML5 shim and Respond.js for IE8 support of HTML5 elements and media queries -->
    <!--[if lt IE 9]>
      <script src="https://oss.maxcdn.com/html5shiv/3.7.2/html5shiv.min.js"></script>
      <script src="https://oss.maxcdn.com/respond/1.4.2/respond.min.js"></script>
    <![endif]-->

	<!-- typekit fonts -->
	<script src="https://use.typekit.net/jtk4slp.js"></script>
	<script>try{Typekit.load({ async: true });}catch(e){}</script>
	<!-- end typekit-->

    <!-- Bootstrap core JavaScript
    ================================================== -->
    <!-- Placed at the end of the document [oops, no the start - jQuery for DataTables] so the pages load faster -->
    <script src="../js/jquery.min.js"></script>
    <script src="../js/bootstrap.min.js"></script>

	<!-- Chart.js -->
	<script src="../Chart.js-master/Chart.js"></script>

	<!-- Highcharts.js -->
	<!-- highcahrts src /Library/WebServer/Documents/Highcharts-4.2.3/js/highcharts.js -->
	<!-- guage -->
	<script src="http://localhost/Highcharts-4.2.3/js/highcharts.js"></script>
	<script src="http://localhost/Highcharts-4.2.3/js/highcharts-more.js"></script>
	<script src="http://localhost/Highcharts-4.2.3/js/modules/solid-gauge.js"></script>
	<!-- heatmap -->
	<!--<script src="https://code.highcharts.com/highcharts.js"></script>-->
	<script src="http://localhost/Highcharts-4.2.3/js/modules/data.js"></script>
	<script src="http://localhost/Highcharts-4.2.3/js/modules/heatmap.js"></script>
	<!-- radar ('polar') -->
	<!--
		<script src="https://code.highcharts.com/highcharts.js"></script>
		<script src="https://code.highcharts.com/highcharts-more.js"></script>
		<script src="https://code.highcharts.com/modules/exporting.js"></script>
	-->
	<script src="http://localhost/Highcharts-4.2.3/js/modules/exporting.js"></script>


	<!-- dataTables -->
	<!--<script type="text/javascript" charset="utf8" src="//cdn.datatables.net/1.10.11/js/jquery.dataTables.js"></script>-->
	<!-- kill the tables - Aug 05 2016
	<script type="text/javascript" charset="utf8" src="../DataTables/DataTables-1.10.11/js/jquery.dataTables.js"></script>
	<script type="text/javascript" charset="utf-8">
		/** KILL THE DATA TABLE RENDER FOR NOW
		$(document).ready(function() {
			$('#blast').dataTable();
		} );
		**/
	</script>
	end kill tables -->
  </head>

  <body>
    <!-- Fixed navbar -->
    <nav class="navbar navbar-default navbar-fixed-top">
      <div class="container">
        <div class="navbar-header">
          <button type="button" class="navbar-toggle collapsed" data-toggle="collapse" data-target="#navbar" aria-expanded="false" aria-controls="navbar">
            <span class="sr-only">Toggle navigation</span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
            <span class="icon-bar"></span>
          </button>
          <a class="navbar-brand" href="#">Real-time phylogenomics</a>
        </div>
        <div id="navbar" class="collapse navbar-collapse">
          <ul class="nav navbar-nav">
            <li class="active"><a href="#">Summary</a></li>
            <li><a href="#BLAST">BLAST</a></li>
            <li><a href="#about">About</a></li>
          </ul>
        </div><!--/.nav-collapse -->
      </div>
    </nav>
    
    <section>
		<div class="container" style="margin-top:50px">
			<!-- Example row of columns -->
			<h1>Summary</h1>
			<div class="row">
				<?php
				//get the basic data file
				echo "<pre style=\"text-size:larger\">Last data: ";
				$mySFfile = fopen("./data/SF.txt", "r") or warn("Unable to open file!");
				while(!feof($mySFfile)) {
					$data = explode("\t",fgets($mySFfile));
						echo "$data[0]\n";
				}
				
				//get the unknown_sample.stats file
				$total_reads = 0;
				$total_yield = 0;
				$mySFfile = fopen("./data/unknown_sample.stats", "r") or warn("Unable to open file!");
				while(!feof($mySFfile)) {
					$data = explode("\t",fgets($mySFfile));
						echo "$data[0] $data[1]";
						if(preg_match('/total\ reads/',$data[0])){
							$total_reads = $data[1];
						}
						if(preg_match('/total\ base\ pairs/',$data[0])){
							$total_yield = $data[1];
						}
				}
				
				//see about those blast hits
				echo "\nBLASTN hits:\n";
				$cum_lengths = array();
				$cum_idents = array();
				$cum_Ns = array();
				$i = 0;
				foreach(array('silene','beta','sorbus','napenthes','erycina','lyrata','thaliana') as $db){
					$mySFfile = fopen("./data/agg.blasts.$db", "r") or warn("Unable to open file!");
					$cum_length = 0;
					$cum_ident = 0;
					$N = 0;
					while(!feof($mySFfile)) {
						$data = explode("\t",fgets($mySFfile));
						if(count($data)>3){
							$cum_length = $data[0];
							$cum_ident = $data[1];
							$N = trim($data[3]);
						}
					}
					echo "| $db: $cum_length $cum_ident $N \n";
					$cum_lengths[$i] = $cum_length;
					$cum_idents[$i] = $cum_ident;
					$cum_Ns[$i] = $N;
					$i++;
				}
				//echo implode(",", $cum_lengths);
				//echo implode(",", $cum_idents);
				//echo implode(",", $cum_Ns);
				echo "</pre>";	
				?>
			</div>
			<div class="row">
				<div class="col-md-4">
					<h2>Read rate (throughput)</h2>
					<p id="container-speed" style="width: 300px; height: 200px"></p>
					<p id="container-rpm" style="width: 300px; height: 200px"></p>
					<p>This shows how many strands of DNA we have read so far, and how many DNA letters we've read overall.</p>
					<p>The number above looks quite big, but we're actually attempting to match DNA with only a <i>tiny</i> fraction of the genome - <b>many plant genomes have more than 1,000,000,000 letters of DNA!</b></p>
					<p><a class="btn btn-default" href="#" role="button">View details &raquo;</a></p>
				</div>
				<div class="col-md-8">
					<h2>ID - Mystery sample 'B'</h2>
					<p id="radar-container"  style="width: 600px"></p>
					<p>This shows which species match the reads we've sequenced. <b>The larger the wedge is, the closer the mystery DNA matches that species.</b></p>
					<p>We're using a process called <b>BLASTN</b>. This is a program that compares the DNA reads we've just sequenced to a collection of reference reads we prepared in the lab last week - a bit like matching pieces in a jigsaw.</p>
					<p><a class="btn btn-default" href="#" role="button">View details &raquo;</a></p>
				</div>
			</div>
		</div>
    </section>
    
    <!-- 
    Flowcell information (may as well write it here):
    
    20160 08 07 AM _ sample_05
    2002012239	FAD22873	1875659	1000 {465,338,157,40}	Sample_D
    2002012435	FAD22885	71007	191 {164,23,4,0}	Sample_B
    
    -->
    <section>
    <a name="BLAST"></a>
    <div class="container">
		<h1>BLAST</h1>
		<div class="row">
			<div class="col-md-4"><h3>Summary table</h3><p>(DataTables render of blast hits)</p>
<!-- DELETED AUG 05 2016 -->
			</div>
		</div>
      </div>
    </section>
    


	<div class="container">
		<footer>
			<p>&copy; Joe Parker 2016</p>
		</footer>
	</div>
	
	
	
	<!-- Highcharts.js RADAR ('Polar') -->
	<script type="text/javascript">
$(function () {

    $('#radar-container').highcharts({

        chart: {
            polar: true
        },

        title: {
            text: 'DNA Matching Engine'
        },

        pane: {
            startAngle: 0,
            endAngle: 360
        },

        xAxis: {
						//categories: ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun'],
						//type: categories,
            tickInterval: 60,
            min: 0,
            max: 360,
             
            labels: {
                
                formatter: function () {
                		var whichTick =  this.value / 60;
                    var myLabels = 
                    	[
                       
                       'Sea Campion',
                       'Beet',
                       'Whitebeam', 
                       'Pitcher Plant',
                       'Dancing Lady Orchid', 
                       'Fuschia'
                       //'Rock-Cress'
                       ];
                    return myLabels[whichTick];
                }
            }
            
            //labels: [1, 2, 3, 4, 5, 6]
        },

        yAxis: {
            min: 0
        },

        plotOptions: {
            series: {
                pointStart: 0,
                pointInterval: 60
            },
            column: {
                pointPadding: 0,
                groupPadding: 0
            }
        },

        series: [{
/*
 *
            type: 'column',
            name: 'Column',
            data: [8, 7, 6, 5, 4, 3, 2],
            pointPlacement: 'between'
        }, {
            type: 'line',
            name: 'Line',
            data: [1, 2, 3, 4, 5, 6, 7, 8]
        }, {
            type: 'area',
            name: 'Area 3',
            data: [1, 12, 3, 10, 1, 7, 4, 5]
        }, {
 *
 */
            type: 'column',
            name: 'Cumulative matched length',
            data: [<?php echo implode(",", $cum_lengths)?>]
        }, {
            type: 'column',
            name: 'Cumulative identities',
            data: [<?php echo implode(",", $cum_idents)?>]
        }]
    });
});	</script>
	<!-- // RADAR -->
	

<!-- Highcharts.js speed guage data -->
<script type=text/javascript>
$(function () {

    var gaugeOptions = {

        chart: {
            type: 'solidgauge'
        },

        title: null,

        pane: {
            center: ['50%', '85%'],
            size: '140%',
            startAngle: -90,
            endAngle: 90,
            background: {
                backgroundColor: (Highcharts.theme && Highcharts.theme.background2) || '#EEE',
                innerRadius: '60%',
                outerRadius: '100%',
                shape: 'arc'
            }
        },

        tooltip: {
            enabled: false
        },

        // the value axis
        yAxis: {
            stops: [
                [0.1, '#55BF3B'], // green
                [0.5, '#DDDF0D'], // yellow
                [0.9, '#DF5353'] // red
            ],
            lineWidth: 0,
            minorTickInterval: null,
            tickAmount: 2,
            title: {
                y: -70
            },
            labels: {
                y: 16
            }
        },

        plotOptions: {
            solidgauge: {
                dataLabels: {
                    y: 5,
                    borderWidth: 0,
                    useHTML: true
                }
            }
        }
    };

    // The speed gauge
    $('#container-speed').highcharts(Highcharts.merge(gaugeOptions, {
        yAxis: {
            min: 0,
            max: 1500,
            title: {
                text: 'Reads'
            }
        },

        credits: {
            enabled: false
        },

        series: [{
            name: 'Speed',
            data: [<?php echo $total_reads ?>],
            dataLabels: {
                format: '<div style="text-align:center"><span style="font-size:25px;color:' +
                    ((Highcharts.theme && Highcharts.theme.contrastTextColor) || 'black') + '">{y}</span><br/>' +
                       '<span style="font-size:12px;color:silver">total reads</span></div>'
            },
            tooltip: {
                valueSuffix: ' km/h'
            }
        }]

    }));

    // The RPM gauge
    $('#container-rpm').highcharts(Highcharts.merge(gaugeOptions, {
        yAxis: {
            min: 0,
            max: 5000000,
            title: {
                text: 'DNA letters'
            }
        },

        series: [{
            name: 'RPM',
            data: [<?php echo $total_yield ?>],
            dataLabels: {
                format: '<div style="text-align:center"><span style="font-size:25px;color:' +
                    ((Highcharts.theme && Highcharts.theme.contrastTextColor) || 'black') + '">{y:.0f}</span><br/>' +
                       '<span style="font-size:12px;color:silver">DNA letters</span></div>'
            },
            tooltip: {
                valueSuffix: ' revolutions/min'
            }
        }]

    }));

    // Bring life to the dials
    setTimeout(function () {
        // Speed
        var chart = $('#container-speed').highcharts(),
            point,
            newVal,
            inc;

        /*
         *
        if (chart) {
            point = chart.series[0].points[0];
            inc = Math.round((Math.random() - 0.5) * 100);
            newVal = point.y + inc;

            if (newVal < 0 || newVal > 200) {
                newVal = point.y - inc;
            }

			//turn off bounciness
            //point.update(newVal);
        }
         *
         */

        // RPM
        /*
         *
        chart = $('#container-rpm').highcharts();
        if (chart) {
            point = chart.series[0].points[0];
            inc = Math.random() - 0.5;
            newVal = point.y + inc;

            if (newVal < 0 || newVal > 5) {
                newVal = point.y - inc;
            }

			//turn off bounciness
            //point.update(newVal);
        }
         *
         */
    }, 2000);


});

</script>
<!-- end guages -->

	<!-- IE10 viewport hack for Surface/desktop Windows 8 bug -->
    <!--<script src="../../assets/js/ie10-viewport-bug-workaround.js"></script>-->
  </body>
</html>
