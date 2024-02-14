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

    return {nodes: new vis.DataView(new vis.DataSet(neededNodes), {filter: filterFunction}),
            edges: new vis.DataSet(neededEdges)};
}


function startNetwork(container, options, predictedCheckboxes, physicsCheckbox,
                      nodeFilterSelector, nodeFilterValue, nodes, edges, graphFilter) {

    initialData = setNodes(nodes, edges, predictedCheckboxes, graphFilter);
    var network2 = new vis.Network(container, initialData, options);

    predictedCheckboxes.forEach(checkbox => {
        checkbox.addEventListener('change', function () {
            // check if length of edges[checkbox.value] is greater than 0, if not do nothing
            if (edges[checkbox.value].length > 0) {
                network2.setData(setNodes(nodes, edges, predictedCheckboxes, graphFilter));
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