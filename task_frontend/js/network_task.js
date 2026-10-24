/**
 * network_task.js
 *
 * The "net" task's client-side logic, extracted verbatim from the live
 * MySocialBrain.org platform's shared js/studylogic.js (which holds the
 * equivalent logic for every study on the platform, not just this one).
 * Shared for transparency about exactly how the task was presented and
 * how choices were recorded -- NOT a standalone, runnable copy. See this
 * folder's README for what's deliberately left out and why.
 *
 * Original location: js/studylogic.js, functions `network` (task trial
 * loop) and `network_results` (post-task feedback visualization), plus
 * the shared helpers `getNetworks`, `renderNetwork`, and `colorBar` that
 * both call into.
 */

// ---------------------------------------------------------------------
// Main task entry point. Called by the platform's page router once the
// consent/instructions flow reaches this study (meta is unused here but
// kept in the signature to match every other study's entry-point shape).
// ---------------------------------------------------------------------
var network = function(meta){
    var EM = window.EM
    var modules = EM.modules;
	var settings = EM.settings;
	var userData = EM.userData;
	var otherData = EM.otherData;
	var functions = EM.functions;
    var tid = window.EM.settings.tid;

    var loadNetwork = function(){
        var response = {'tid': tid};
        var nnet = 30;
        getNetworks(nnet).done(function(response){
            var networks = response;
            modules.getPage(tid + "_test.html", "modal-body").done(function(response){
                modules.toggleClass("progarea", "hidden");
                runNetwork(networks);
            });
        });
    };

    var runNetwork = function(networks){

        // initialize variables
        var d = new Date();
        var start = d.getTime();
        var curq = 1;
        var results = [];
        var centrality = {"degree":[],"closeness":[],"betweenness":[],"eigenvector":[]};
        var cmean = {"degree":4.29,"closeness":.0196,"betweenness":4.38,"eigenvector":.0574};
        var csd = {"degree":1.68,"closeness":.00381,"betweenness":21.97,"eigenvector":.3456};
        var prog = 0;
        var qn = networks.length;
        var selectedNode;

        // render first network
        $("#submit").attr('disabled',true);
        var unode = modules.randInt(0,networks[0].nodes.length-1);
        renderNetwork(networks[curq-1],unode,'test');
        $(".node").on('click tap',function(e){
           e.preventDefault();
           if (e.target.classList.length == 1){
               d3.selectAll(".selectednode").classed("selectednode",false);
               selectedNode = parseInt(e.target.id);
               var circles = d3.selectAll("circle");
               circles[0][selectedNode].classList.add("selectednode");
               $("#submit").attr('disabled',false);
           };
        });

        // begin trial sequence
        $("#submit").click(function(e){
                var d = new Date();
                var rt = d.getTime() - start;
                var prediction = response.data;
                Object.keys(centrality).forEach(function(ctype){
                    centrality[ctype].push(networks[curq-1].nodes[selectedNode][ctype]);
                });
                results.push({'netid':networks[curq-1].id,'choiceid':selectedNode,'rt':rt,'unode':unode});
                if(curq == qn){
                    $("#submit").attr('disabled',true);
                    if(userData.userID){
                        var id = userData.userID;
                    } else {
                        var id = userData.guestID;
                    }
                    var resdat = {"tid":tid,"id":id,"data":results};
                    modules.post(resdat,"netresponses.php",null);
                    Object.keys(centrality).forEach(function(ctype){
                        centrality[ctype] = (modules.mean(centrality[ctype]) - cmean[ctype])/csd[ctype];
                    });
                    otherData.visdata = centrality;
                    functions.initPost();
                    $('#progress-bar').width(0 + "%");
                    modules.toggleClass("progarea", "hidden");
                } else {
                    $("#submit").attr('disabled',true);
                    prog = Math.round(100*(curq)/qn);
                    $('#progress-bar').width(prog + "%");
                    $("#visarea").empty();
                    d = new Date();
                    start = d.getTime();
                    curq = ++curq;
                    unode = modules.randInt(0,networks[curq-1].nodes.length-1);
                    renderNetwork(networks[curq-1],unode,'test');
                    $(".node").on('click tap',function(e){
                        e.preventDefault();
                       if (e.target.classList.length == 1){
                           d3.selectAll(".selectednode").classed("selectednode",false);
                           selectedNode = parseInt(e.target.id);
                           var circles = d3.selectAll("circle");
                           circles[0][selectedNode].classList.add("selectednode");
                           $("#submit").attr('disabled',false);
                       };
                    });
                }
        });
    };

    loadNetwork();
};


// ---------------------------------------------------------------------
// Post-task feedback: per-centrality-type "you vs. average" bar chart,
// plus an interactive example network participants can re-color by
// centrality type.
// ---------------------------------------------------------------------
var network_results = function(){
    var EM = window.EM
    var modules = EM.modules;
	var settings = EM.settings;
	var userData = EM.userData;
	var otherData = EM.otherData;
	var functions = EM.functions;
    tid = window.EM.settings.tid;

    modules.getPage(tid + '_results.html','modal-body').done(function(response){
        functions.initShare();
        if(document.getElementById("debrief_general") != null){
            $('#debrief_general').load("index.html #rcontainer");
        }
        if(document.getElementById("modalcontactinsert") != null){
            $('#modalcontactinsert').load("pages/contact.html #contact");
        }
        modules.post(response,"getvis.php",this.id).complete(function(response){
            var mcen = JSON.parse(response.responseJSON["data"][0]["visdata"]);
            if (otherData.visdata == null){
                response = null;
                modules.post(response,"getindivis.php",this.id).complete(function(response){
                    response = JSON.parse(response.responseJSON.data);
                    var centrality = {};
                    Object.keys(response).forEach(function(key) {
                        centrality[key] = parseFloat(response[key]);
                    });
                    makeBars(centrality,mcen);
                });
            } else {;
                makeBars(otherData.visdata,mcen);
            }
        });
    });

    var makeBars = function(centrality,mcen){

        // organize data
        var data = [];
        var ind = 0;
        var allvals = []
        Object.keys(centrality).forEach(function(key) {
            data.push({dim:key,val:centrality[key],type:" (you)",ind:ind});
            allvals.push(centrality[key]);
            ind+=1
            data.push({dim:key,val:mcen[key]/mcen["count"],type:" (average)",ind:ind});
            allvals.push(mcen[key]/mcen["count"]);
            ind+=1
        });
        var nd = data.length/2;

        // plot margins
        var margin = {top: 60, right: 50, bottom: 10, left: 120},
            width = Math.min(500,window.innerWidth/2),
            height = nd*100;

        // graph outline
        var svg = d3.select("#bararea").append("svg")
            .attr("width", width + margin.left + margin.right)
            .attr("height", height + margin.top + margin.bottom)
          .append("g")
            .attr("transform", "translate(" + margin.left + "," + margin.top + ")");


        // scales
        var x = d3.scale.linear()
            .domain([Math.min(d3.min(allvals),0),d3.max(allvals)])
            .range([0,width]);

        var y = d3.scale.ordinal()
            .domain(data.map(function(d) { return d.dim + d.type; }))
            .rangeRoundBands([0, height], .2);

        var color = d3.scale.category20();

        // create bars
        svg.selectAll(".bar")
            .data(data)
          .enter().append("rect")
            .attr("class", "bar")
            .attr("x", function(d) { return x(Math.min(0, d.val)); })
            .attr("y", function(d) { return y(d.dim + d.type); })
            .attr("width", function(d) { return Math.abs(x(d.val) - x(0)); })
            .attr("height", y.rangeBand())
            .attr("fill", function(d) { return color(d.ind); } );

        // axes
        var xAxis = d3.svg.axis()
            .scale(x)
            .ticks(5)
            .orient("top");

        var yAxis = d3.svg.axis()
            .scale(y)
            .tickPadding(x(-.9))
            .orient("left");

        svg.append("g")
            .attr("class", "y axis")
            .attr("transform", "translate(" + x(0) + ",0)")
            .call(yAxis);

         svg.append("g")
            .attr("class", "x axis")
            .call(xAxis);

        svg.append("text")
            .attr("class", "xlabel")
            .attr("text-anchor", "middle")
            .attr("x", width/2)
            .attr("y", -30)
            .text("Weight given to each type of centrality:");

        getNetworks(1).done(function(response){
            var networks = response;
            renderNetwork(networks[0],0,"results");
        });

    };
};


// ---------------------------------------------------------------------
// Shared helpers used by both network() and network_results() above.
// ---------------------------------------------------------------------

// Fetches `nnet` random stimulus graphs from the server-side pool (see
// generate_networks.R in the main analysis repository for how that pool
// of 10,000 graphs was generated) by random integer ID; the backend
// endpoint get_networks.php is not included here (server-side).
var getNetworks = function(nnet){
    var deferred = $.Deferred();
    var modules = window.EM.modules;
    var response = {"nnet": nnet};
    var netids = [];
    for (var i = 1; i <= nnet; i ++){
        netids.push(modules.randInt(1,10000));
    }
    response = {"netids": netids};
    modules.post(response, "get_networks.php", null).complete(function(response){
        response = response.responseJSON["data"];
        var networks = [];
        for (var i = 0; i < nnet; i ++){
            networks.push(JSON.parse(response[i][0]["graph"]));
        }
        deferred.resolve(networks);
    });
    return deferred;
};

// D3 force-directed rendering of one network into #visarea. vistype
// 'test' marks the participant's own node and existing friends and wires
// up node selection; vistype 'results' instead adds the centrality
// color-bar legend and lets the viewer toggle which centrality measure
// colors the nodes.
var renderNetwork = function(network,usernode,vistype){
    // vistypes = test, results

    var visarea = d3.select("#visarea");
    var width = visarea.node().getBoundingClientRect().width;
    var height = Math.max(window.innerHeight*.5,300);
    var radius = Math.max(height/40,10);
    var charge = -525*(Math.log10(width)/4);

    var svg = d3.select("#visarea").append("svg")
        .attr("width", width)
        .attr("height", height);

    var force = d3.layout.force()
        .nodes(network.nodes)
        .links(network.links)
        .size([width, height])
        .charge(charge)
        .alpha(.05)
        .on("tick", tick)
        .start();

    var link = svg.selectAll(".link")
       .data(network.links)
     .enter().append("line")
       .attr("class", "link");

    var node = svg.selectAll(".node")
       .data(network.nodes)
     .enter().append("circle")
       .attr("id",function(d) { return d.id; })
       .attr("class", "node")
       .attr("r", radius);

    function tick() {

      node.attr("cx", function(d) { return Math.max(Math.min(d.x,width-radius),radius); })
          .attr("cy", function(d) { return Math.max(Math.min(d.y,height-radius),radius); });

      link.attr("x1", function(d) { return Math.max(Math.min(d.source.x,width-radius),radius); })
          .attr("y1", function(d) { return Math.max(Math.min(d.source.y,height-radius),radius); })
          .attr("x2", function(d) { return Math.max(Math.min(d.target.x,width-radius),radius); })
          .attr("y2", function(d) { return Math.max(Math.min(d.target.y,height-radius),radius); });

    }

    var circles = d3.selectAll("circle");

    // map color to network
    var mapColor = function(ctype){
        var censcale = d3.scale.linear()
            .domain([d3.min(network.nodes.map(function(d){return d[ctype];})), d3.max(network.nodes.map(function(d){return d[ctype];}))])
            .range(["black", "red"]);
        circles.style("fill",function(d) {return censcale(d[ctype]);});
    }

    // change behavior depending on vistype
    if (vistype == "test"){
        circles[0][usernode].classList.add("usernode");
        var nlink = network.links.length;
        for (var i = 0; i < nlink; i ++){
            if (network.links[i].source.id == usernode){
                circles[0][network.links[i].target.id].classList.add("friendnode");
            } else if (network.links[i].target.id == usernode) {
                circles[0][network.links[i].source.id].classList.add("friendnode");
            }
        }
    } else {
        // add color bar
        $("#visarea").append("<div id='colorbar'></div>");
        colorBar("#colorbar",width*.6,40,"black","red");

        // map initial centrality
        mapColor("degree");

        // monitor and update centrality mapping
        $(".ctoggle").click(function(e){
             d3.selectAll(".ctoggle").classed("active",false);
             sel = e.target.id;
             $("#" + sel).addClass("active");
             mapColor(sel);
        });


    }

}

// Simple linear-gradient legend bar for the results-page centrality
// color scale.
var colorBar = function(target,width,height,startcol,endcol){

    var svg = d3.select(target).append('svg:svg')
        .attr('width',width)
        .attr('height',height);

    var gradient = svg.append("svg:defs")
            .append("svg:linearGradient")
                        .attr("id", "gradient")
                        .attr("x1", "0%")
                        .attr("y1", "0%")
                        .attr("x2", "100%")
                        .attr("y2", "0%")
                        .attr("spreadMethod", "pad");

                    gradient.append("svg:stop")
                        .attr("offset", "0%")
                        .attr("stop-color", startcol)
                        .attr("stop-opacity", 1);

                    gradient.append("svg:stop")
                        .attr("offset", "100%")
                        .attr("stop-color", endcol)
                        .attr("stop-opacity", 1);

                    svg.append("svg:rect")
                        .attr("width", width)
                        .attr("height", height)
                        .style("fill", "url(#gradient)");
};
