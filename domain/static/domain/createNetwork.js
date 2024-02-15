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


function startNetwork(container, options, predictedCheckboxes, physicsCheckbox,
                      nodeFilterSelector, nodeFilterValue, nodes, edges, graphFilter) {

    const initialData = setNodes(nodes, edges, predictedCheckboxes, graphFilter);
    let previousData = initialData;
    var network2 = new vis.Network(container, initialData, options);

    predictedCheckboxes.forEach(checkbox => {
        checkbox.addEventListener('change', function () {
            // check if length of edges[checkbox.value] is greater than 0, if not do nothing
            if (edges[checkbox.value].length > 0) {
                data = setNodes(nodes, edges, predictedCheckboxes, graphFilter);
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
    network2.setData(setNodes(nodes, edges, predictedCheckboxes, graphFilter));
    });
    network2.fit();
}