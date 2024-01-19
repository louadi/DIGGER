function addOptionsToSelect(data) {
    // Get all elements that start with org_select
    var selectElements = document.querySelectorAll('[id^="org_select"]');

    // Add options from 'organisms' and use 'trivial_names' as values
    data.organisms.forEach(function (organism, index) {
        var optionElement = document.createElement('option');
        optionElement.value = data.trivial_names[index];
        optionElement.textContent = organism;

        // Add options to all select elements
        selectElements.forEach(function (selectElement) {
            var clonedOptionElement = optionElement.cloneNode(true);
            if (data.trivial_names[index] === 'human') {
                clonedOptionElement.selected = true;
            }
            selectElement.appendChild(clonedOptionElement);
        });
    });
}

document.addEventListener('DOMContentLoaded', function () {
    fetch('/get_organisms/')
        .then(response => response.json())
        .then(data => {
            addOptionsToSelect(data);
        })
        .catch(error => console.error('Error fetching organisms:', error));
});