// Description: This file contains the functions to create the network visualizations using vis.js

// Set the vis.js dataset based on the predictedCheckboxes and the filterFunction
function setNodes(nodes, edges, predictedCheckboxes, filterFunction) {
    neededNodes = nodes['original'].slice()
    neededEdges = edges['original'].slice()

    let seenIds = new Set();
    nodes['original'].forEach(node => {
        seenIds.add(node.id);
    });

    predictedCheckboxes.forEach(checkbox => {
        if (checkbox.checked) {
            nodes[checkbox.value].forEach(node => {
                if (!seenIds.has(node.id)) {
                    neededNodes.push(node);
                    seenIds.add(node.id);
                }
            });
            neededEdges = neededEdges.concat(edges[checkbox.value]);
        }
    })
    // change this to filter neededNodes by the filterFunction because DataView sucks
    nodes = new vis.DataSet(neededNodes);
    nodes = new vis.DataSet(nodes.get({filter: filterFunction}));

    return {nodes: nodes,
            edges: new vis.DataSet(neededEdges)};
}


const createSetFromEdges = (data, edges) => {
    const edgeSet = new Set();
    edges.forEach(id => {
        const edge = data.edges.get(id);
        // skip if edge is not in nodes ids
        if (!data.nodes.get(edge.from) || !data.nodes.get(edge.to)) {
            return;
        }
        edgeSet.add(JSON.stringify([edge.from, edge.to]));
    });
    return edgeSet;
};


function compareDataSets(data1, data2) {
    // immediately return false if the length of the nodes are different
    if (data1.nodes.length !== data2.nodes.length) {
        return false;
    }
    // check if data1.nodes ids are in data2.nodes ids
    data1.nodes.getIds().forEach(id => {
        if (!data2.nodes.get(id)) {
            return false;
        }
    });
    // add from and to of each edge as a tuple to a set
    const data1_edges = createSetFromEdges(data1, data1.edges.getIds());
    const data2_edges = createSetFromEdges(data2, data2.edges.getIds());
    let intersect = new Set([...data1_edges].filter(i => data2_edges.has(i)));
    // check if the set of edges are the same
    return intersect.size === data1_edges.size && intersect.size === data2_edges.size;

}

// unfortunately, ProteinView must be filtered in the frontend to avoid floating nodes while retaining floating nodes
// in the InteractionView. This is a trade-off to keep the sent data small
function filterProteinView(networkData, filterValue) {
    if (filterValue !== "PPI" && filterValue !== "PPI-DDI") {
        return networkData;
    }
    nodes = networkData.nodes.get();
    edges = networkData.edges.get();

    if (filterValue === "PPI-DDI") {
        // filter nodes that contain / in their id,
        depthNotReached = nodes.filter(node => connectedToMainNetwork(edges, node.id));
        nodes = nodes.filter(node => !depthNotReached.includes(node));
    } else {
        // get all the ids of the nodes that are in edges
        let seenIds = new Set();
        edges.forEach(edge => {
            // skip edges with / in them
            if (edge.from.includes("/") && filterValue === 'PPI') {
                return;
            }
            seenIds.add(edge.from);

            if (edge.to.includes("/") && filterValue === 'PPI') {
                return;
            }
            seenIds.add(edge.to)

        });
        // filter out nodes that are not in edges
        nodes = nodes.filter(node => seenIds.has(node.id));
    }

    return {nodes: new vis.DataSet(nodes), edges: networkData.edges};
}


// function that checks if a node is connected to the main network (to avoid random floating nodes)
function connectedToMainNetwork(edges, startNode) {
    const visited = new Set();
    let depthReached = false;
    // I'm not proud of what I'm about to do
    // For nodes that are a Domain belonging to a protein, the depth is 2 since proteins with >1 domain allow for a
    // depth of 2, however if the protein is connected to the main network, it must at least be 3
    const requiredDepth = startNode.includes("/") ? 2 : 1;


    function dfs(node, depth) {
        visited.add(node);
        if (depth > requiredDepth) {
            depthReached = true;
            return true;
        }
        edges.forEach(edge => {
            if (edge.from === node && !visited.has(edge.to)) {
                if (dfs(edge.to, depth + 1)) {
                    return true;
                }
            }
            if (edge.to === node && !visited.has(edge.from)) {
                if (dfs(edge.from, depth + 1)) {
                    return true;
                }
            }
        });
        return false;
    }

    dfs(startNode, 0);
    return !depthReached;
}


// start the network as dynamically as possible to avoid code duplication
function startNetwork(container, options, predictedCheckboxes, physicsCheckbox,
                      nodeFilterSelector, nodeFilterValue, nodes, edges, graphFilter) {

    let initialData = setNodes(nodes, edges, predictedCheckboxes, graphFilter);
    if (nodeFilterValue.value === "PPI" || nodeFilterValue.value === "PPI-DDI") {
        initialData = filterProteinView(initialData, nodeFilterValue.value);
    }
    let previousData = initialData;
    var network2 = new vis.Network(container, initialData, options);

    predictedCheckboxes.forEach(checkbox => {
        checkbox.addEventListener('change', function () {
            // check if length of edges[checkbox.value] is greater than 0, if not do nothing
            if (edges[checkbox.value].length > 0) {
                data = setNodes(nodes, edges, predictedCheckboxes, graphFilter);
                if (nodeFilterValue.value === "PPI" || nodeFilterValue.value === "PPI-DDI") {
                    data = filterProteinView(data, nodeFilterValue.value);
                }
                // check if the nodes and edges are the same as the previous ones, if not update the network
                if (!compareDataSets(data, previousData)) {
                    previousData = data;
                    network2.setData(data);
                }
            }
        });
    });

    physicsCheckbox.addEventListener('change', function () {
        if (this.checked) {
            network2.setOptions({physics: true});
        } else {
            network2.setOptions({physics: false});
        }
    });

    nodeFilterSelector.addEventListener("change", function(e) {
    // set new value to filter variable
    nodeFilterValue.value = e.target.value;
    data = setNodes(nodes, edges, predictedCheckboxes, graphFilter);
    if (nodeFilterValue.value === "PPI" || nodeFilterValue.value === "PPI-DDI") {
        data = filterProteinView(data, nodeFilterValue.value);
    }
    network2.setData(data);
    });
    network2.fit();
}